#include "Coincidences.h"

// Implementation of the low-level coincidence counting logic. Keeping detailed
// comments here helps both the CLI driver and the Python wrapper stay in sync
// about expectations (picosecond arithmetic, span-based views, etc.).

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
struct DelayScanConfig {
  // All members are stored in picoseconds to avoid repeated conversions.
  long long startPs = 0;
  long long endPs = 0;
  long long stepPs = 0;
  size_t steps = 0;
};

DelayScanConfig buildConfig(long long delayStartPs, long long delayEndPs,
                            long long delayStepPs) {
  if (delayStepPs <= 0)
    throw std::invalid_argument("delayStep must be positive in ps");

  if (delayEndPs < delayStartPs)
    return {delayStartPs, delayEndPs, delayStepPs, 0};

  const size_t steps =
      static_cast<size_t>(((delayEndPs - delayStartPs) / delayStepPs) + 1);
  return {delayStartPs, delayEndPs, delayStepPs, steps};
}
} // namespace

std::span<const long long>
appendNextFirstEvent(std::span<const long long> currentSecond,
                     std::span<const long long> nextSecond,
                     std::vector<long long> &scratch) {
  // Fast-path: when there's nothing in the next bucket we can just return a
  // span over the original memory—no copies, no allocations.
  if (nextSecond.empty()) {
    if (currentSecond.empty()) {
      scratch.clear();
      return {};
    }
    return currentSecond;
  }

  // Slow-path: need to append the head of the next bucket to preserve
  // possible coincidences that straddle the second boundary.
  scratch.assign(currentSecond.begin(), currentSecond.end());
  scratch.push_back(nextSecond.front());
  return std::span<const long long>(scratch.data(), scratch.size());
}

namespace {
template <bool Collect>
int countCoincidencesWithDelay(
    std::span<const long long> ch1, std::span<const long long> ch2,
    long long coincWindowPs, long long delayPs,
    std::vector<std::pair<long long, long long>> *outHits) {
  const long long halfWindow = coincWindowPs / 2;
  const long long lowerBound = -halfWindow;
  const long long upperBound = halfWindow;

  int count = 0;
  size_t i = 0;
  size_t j = 0;
  const size_t size1 = ch1.size();
  const size_t size2 = ch2.size();

  while (i < size1 && j < size2) {
    const long long shifted = ch1[i] - delayPs;
    const long long diff = shifted - ch2[j];

    if (diff < lowerBound) {
      ++i;
    } else if (diff > upperBound) {
      ++j;
    } else {
      ++count;
      if constexpr (Collect) {
        outHits->emplace_back(ch1[i], ch2[j]);
      }
      ++i;
      ++j;
    }
  }
  return count;
}
} // namespace

std::vector<std::pair<long long, long long>>
collectCoincidencesWithDelay(std::span<const long long> ch1,
                             std::span<const long long> ch2,
                             long long coincWindowPs, long long delayPs) {
  std::vector<std::pair<long long, long long>> hits;
  countCoincidencesWithDelay<true>(ch1, ch2, coincWindowPs, delayPs, &hits);
  return hits;
}

int countCoincidencesWithDelay(std::span<const long long> ch1,
                               std::span<const long long> ch2,
                               long long coincWindowPs, long long delayPs) {
  return countCoincidencesWithDelay<false>(ch1, ch2, coincWindowPs, delayPs,
                                           nullptr);
}

int countNFoldCoincidences(
    const std::vector<std::span<const long long>> &channels,
    long long coincWindowPs, std::span<const long long> offsetsPs) {
  if (channels.size() < 2)
    throw std::invalid_argument(
        "At least two channels required for coincidences");
  if (!offsetsPs.empty() && offsetsPs.size() != channels.size())
    throw std::invalid_argument("offsets size must match channels size");
  if (channels.size() == 2 && offsetsPs.empty())
    return countCoincidencesWithDelay(channels[0], channels[1], coincWindowPs,
                                      0);

  struct Tagged {
    long long timestamp;
    size_t channelIdx;
  };

  size_t totalEvents = 0;
  for (auto span : channels)
    totalEvents += span.size();
  if (totalEvents == 0)
    return 0;

  std::vector<Tagged> merged;
  merged.reserve(totalEvents);
  for (size_t idx = 0; idx < channels.size(); ++idx) {
    const long long offset = offsetsPs.empty() ? 0 : offsetsPs[idx];
    for (long long ts : channels[idx])
      merged.push_back({ts + offset, idx});
  }
  std::sort(merged.begin(), merged.end(), [](const Tagged &a, const Tagged &b) {
    return a.timestamp < b.timestamp;
  });

  std::vector<int> freq(channels.size(), 0);
  size_t have = 0;
  size_t left = 0;
  int coincidences = 0;

  for (size_t right = 0; right < merged.size(); ++right) {
    const size_t idx = merged[right].channelIdx;
    if (++freq[idx] == 1)
      ++have;

    while (merged[right].timestamp - merged[left].timestamp > coincWindowPs &&
           left < right) {
      const size_t lidx = merged[left].channelIdx;
      if (--freq[lidx] == 0)
        --have;
      ++left;
    }

    if (have == channels.size()) {
      ++coincidences;
      const size_t lidx = merged[left].channelIdx;
      if (--freq[lidx] == 0)
        --have;
      ++left;
    }
  }

  return coincidences;
}

DelayScan scanDelays(std::span<const long long> channel1,
                     std::span<const long long> channel2,
                     long long coincWindowPs,
                     long long delayStartPs, long long delayEndPs,
                     long long delayStepPs) {
  DelayScan scan;
  const DelayScanConfig config =
      buildConfig(delayStartPs, delayEndPs, delayStepPs);
  if (config.steps == 0)
    return scan;

  scan.delaysPs.resize(config.steps);
  scan.counts.assign(config.steps, 0);
  for (size_t idx = 0; idx < config.steps; ++idx)
    scan.delaysPs[idx] = static_cast<double>(
        config.startPs + static_cast<long long>(idx) * config.stepPs);

  if (channel1.empty() || channel2.empty())
    return scan; // counts already zero-initialized

  // Difference array (size = steps + 1 so "end + 1" stays in-bounds).
  std::vector<long long> diff(config.steps + 1, 0);
  size_t jLo = 0;
  size_t jHi = 0;
  const long long halfWindow = coincWindowPs / 2;
  const long long minNeeded = config.startPs - halfWindow;
  const long long maxNeeded = config.endPs + halfWindow;

  for (const long long t1 : channel1) {
    // Keep channel2[jLo:jHi) aligned with timestamps that can still
    // contribute coincidences for this t1 once the delay range is applied.
    const long long lowCut = t1 - maxNeeded;
    while (jLo < channel2.size() && channel2[jLo] < lowCut)
      ++jLo;

    const long long highCut = t1 - minNeeded;
    if (jHi < jLo)
      jHi = jLo;
    while (jHi < channel2.size() && channel2[jHi] <= highCut)
      ++jHi;

    for (size_t j = jLo; j < jHi; ++j) {
      const long long diffCenter = t1 - channel2[j];
      long long intervalStart = diffCenter - halfWindow;
      long long intervalEnd = diffCenter + halfWindow;
      if (intervalEnd < config.startPs || intervalStart > config.endPs)
        continue;
      intervalStart = std::max(intervalStart, config.startPs);
      intervalEnd = std::min(intervalEnd, config.endPs);

      // Round into discrete delay bins: start indexes the first bin whose
      // centre lies inside the window; end indexes the last bin.
      const long long offsetStart = intervalStart - config.startPs;
      const long long offsetEnd = intervalEnd - config.startPs;
      const size_t idxStart = static_cast<size_t>(
          (offsetStart + config.stepPs - 1) / config.stepPs);
      const size_t idxEnd = static_cast<size_t>(offsetEnd / config.stepPs);
      if (idxStart > idxEnd || idxEnd >= config.steps)
        continue;

      diff[idxStart] += 1;
      diff[idxEnd + 1] -= 1;
    }
  }

  // Prefix-sum the diff array to convert it into actual coincidence counts.
  long long running = 0;
  for (size_t idx = 0; idx < config.steps; ++idx) {
    running += diff[idx];
    scan.counts[idx] = running;
  }
  return scan;
}

BestDelayResult findBestDelay(std::span<const long long> reference,
                              std::span<const long long> target,
                              long long coincWindowPs,
                              long long delayStartPs, long long delayEndPs,
                              long long delayStepPs) {
  BestDelayResult result;
  result.scan = scanDelays(reference, target, coincWindowPs, delayStartPs,
                           delayEndPs, delayStepPs);
  result.bestDelayPs = delayStartPs;
  long long bestCount = std::numeric_limits<long long>::min();
  for (size_t idx = 0; idx < result.scan.counts.size(); ++idx) {
    if (result.scan.counts[idx] > bestCount) {
      bestCount = result.scan.counts[idx];
      result.bestDelayPs =
          static_cast<long long>(std::llround(result.scan.delaysPs[idx]));
    }
  }
  return result;
}

namespace {
// Picks the delay at the centre of the widest run of bins tied for the
// highest count. Returns nullopt when the scan is empty or every bin is zero.
std::optional<double> bestDelayFromScan(const DelayScan &scan) {
  if (scan.delaysPs.empty())
    return std::nullopt;

  const long long best = *std::max_element(scan.counts.begin(), scan.counts.end());
  if (best <= 0)
    return std::nullopt;

  size_t bestRunStart = 0;
  size_t bestRunLen = 0;
  size_t i = 0;
  while (i < scan.counts.size()) {
    if (scan.counts[i] != best) {
      ++i;
      continue;
    }
    const size_t start = i;
    while (i < scan.counts.size() && scan.counts[i] == best)
      ++i;
    if (i - start > bestRunLen) {
      bestRunStart = start;
      bestRunLen = i - start;
    }
  }

  const size_t lastIdx = bestRunStart + bestRunLen - 1;
  return (scan.delaysPs[bestRunStart] + scan.delaysPs[lastIdx]) / 2.0;
}
} // namespace

TwoStageDelayResult findBestDelayTwoStage(std::span<const long long> reference,
                                          std::span<const long long> target,
                                          long long windowPs,
                                          long long coarseWindowPs,
                                          long long coarseHalfRangePs,
                                          long long coarseStepPs,
                                          long long fineHalfRangePs,
                                          long long fineStepPs) {
  TwoStageDelayResult result;

  const DelayScan coarseScan =
      scanDelays(reference, target, coarseWindowPs, -coarseHalfRangePs,
                coarseHalfRangePs, coarseStepPs);
  const std::optional<double> coarseBest = bestDelayFromScan(coarseScan);
  if (!coarseBest)
    return result;

  const long long center = static_cast<long long>(std::llround(*coarseBest));
  result.fineScan = scanDelays(reference, target, windowPs,
                               center - fineHalfRangePs,
                               center + fineHalfRangePs, fineStepPs);
  const std::optional<double> fineBest = bestDelayFromScan(result.fineScan);
  result.bestDelayPs = fineBest ? fineBest : coarseBest;
  return result;
}

WindowEstimate estimateCoincidenceWindow(std::span<const long long> channel1,
                                         std::span<const long long> channel2,
                                         long long delayPs, long long spanPs,
                                         long long binPs) {
  WindowEstimate result;
  // A window of one bin width means each pair lands in essentially the one
  // bin closest to its true (channel1 - channel2) difference, giving a
  // near-unsmoothed histogram of the raw timing peak rather than the
  // coincidence *counts at a window* that scanDelays is normally used for.
  // Clamp so the half-width scanDelays derives from this is never zero.
  const long long binWindowPs = std::max<long long>(2, binPs);
  result.histogram =
      scanDelays(channel1, channel2, binWindowPs, delayPs - spanPs,
                delayPs + spanPs, binPs);
  const std::vector<long long> &counts = result.histogram.counts;
  if (counts.empty())
    return result;

  const size_t peakIdx = static_cast<size_t>(
      std::max_element(counts.begin(), counts.end()) - counts.begin());
  result.peakCount = counts[peakIdx];

  // Background: median count in the outer eighth of bins on each side, as a
  // robust stand-in for the accidental-coincidence floor.
  const size_t edge = std::max<size_t>(1, counts.size() / 8);
  std::vector<long long> tails;
  tails.insert(tails.end(), counts.begin(), counts.begin() + static_cast<long>(edge));
  tails.insert(tails.end(), counts.end() - static_cast<long>(edge), counts.end());
  std::sort(tails.begin(), tails.end());
  result.backgroundCount = static_cast<double>(tails[tails.size() / 2]);

  if (result.peakCount <= result.backgroundCount)
    return result; // no discernible peak above background

  const double halfMax =
      result.backgroundCount +
      (static_cast<double>(result.peakCount) - result.backgroundCount) / 2.0;

  size_t lo = peakIdx;
  while (lo > 0 && static_cast<double>(counts[lo]) >= halfMax)
    --lo;
  size_t hi = peakIdx;
  while (hi + 1 < counts.size() && static_cast<double>(counts[hi]) >= halfMax)
    ++hi;

  result.fwhmPs = result.histogram.delaysPs[hi] - result.histogram.delaysPs[lo];
  result.windowPs = result.fwhmPs; // window is a full width, same as fwhmPs
  return result;
}

void writeResultsToFile(const DelayScan &scan, const std::string &filename) {
  std::ofstream out(filename);
  if (!out.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }
  for (size_t idx = 0; idx < scan.delaysPs.size(); ++idx)
    out << scan.delaysPs[idx] << "," << scan.counts[idx] << "\n";
}
