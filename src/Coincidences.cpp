#include "Coincidences.h"

// Implementation of the low-level coincidence counting logic. Keeping detailed
// comments here helps both the CLI driver and the Python wrapper stay in sync
// about expectations (picosecond arithmetic, span-based views, etc.).

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
constexpr long long kPicosecondsPerNanosecond = 1000LL;

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
appendNextFirstEvent(const std::vector<long long> &currentSecond,
                     const std::vector<long long> &nextSecond,
                     std::vector<long long> &scratch) {
    // Fast-path: when there's nothing in the next bucket we can just return a
    // span over the original memory—no copies, no allocations.
    if (nextSecond.empty()) {
        if (currentSecond.empty()) {
            scratch.clear();
            return {};
        }
        return std::span<const long long>(currentSecond.data(),
                                          currentSecond.size());
    }

    // Slow-path: we need to append the head of the next bucket to preserve
    // possible coincidences that straddle the second boundary.
    scratch.assign(currentSecond.begin(), currentSecond.end());
    scratch.push_back(nextSecond.front());
    return std::span<const long long>(scratch.data(), scratch.size());
}

int countCoincidencesWithDelay(std::span<const long long> ch1,
                               std::span<const long long> ch2,
                               long long coincWindowPs, long long delayPs) {
    // The window is symmetric around zero and expressed as a half-width to keep
    // the math tight.
    const long long lowerBound = -coincWindowPs;
    const long long upperBound = coincWindowPs;

    int count = 0;
    size_t i = 0;
    size_t j = 0;
    const size_t size1 = ch1.size();
    const size_t size2 = ch2.size();

    while (i < size1 && j < size2) {
        const long long shifted = ch1[i] - delayPs;
        const long long diff = shifted - ch2[j];

        // Classic two-pointer sweep over sorted timestamps.
        if (diff < lowerBound) {
            ++i;
        } else if (diff > upperBound) {
            ++j;
        } else {
            ++count;
            ++i;
            ++j;
        }
    }
    return count;
}

int countNFoldCoincidences(const std::vector<std::span<const long long>> &channels,
                           long long coincWindowPs,
                           std::span<const long long> offsetsPs) {
    if (channels.size() < 2)
        throw std::invalid_argument("At least two channels required for coincidences");
    if (!offsetsPs.empty() && offsetsPs.size() != channels.size())
        throw std::invalid_argument("offsets size must match channels size");
    if (channels.size() == 2 && offsetsPs.empty())
        return countCoincidencesWithDelay(channels[0], channels[1], coincWindowPs, 0);

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
    std::sort(merged.begin(), merged.end(),
              [](const Tagged &a, const Tagged &b) { return a.timestamp < b.timestamp; });

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

void computeCoincidencesForRange(std::span<const long long> channel1,
                                 std::span<const long long> channel2,
                                 long long coincWindowPs,
                                 long long delayStartPs, long long delayEndPs,
                                 long long delayStepPs,
                                 std::vector<std::pair<float, int>> &results) {
    results.clear();
    const DelayScanConfig config =
        buildConfig(delayStartPs, delayEndPs, delayStepPs);
    if (config.steps == 0)
        return;

    results.resize(config.steps);
    if (channel1.empty() || channel2.empty()) {
        // Exit when either channel has no events.
        for (size_t idx = 0; idx < config.steps; ++idx) {
            const long long delayPs =
                config.startPs + static_cast<long long>(idx) * config.stepPs;
            results[idx] = {static_cast<float>(delayPs) /
                                kPicosecondsPerNanosecond,
                            0};
        }
        return;
    }

    // Difference array (size = steps + 1 so "end + 1" stays in-bounds).
    std::vector<long long> diff(config.steps + 1, 0);
    size_t jLo = 0;
    size_t jHi = 0;
    const long long minNeeded = config.startPs - coincWindowPs;
    const long long maxNeeded = config.endPs + coincWindowPs;

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
            long long intervalStart = diffCenter - coincWindowPs;
            long long intervalEnd = diffCenter + coincWindowPs;
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
        const long long delayPs =
            config.startPs + static_cast<long long>(idx) * config.stepPs;
        results[idx] = {static_cast<float>(delayPs) /
                            kPicosecondsPerNanosecond,
                        static_cast<int>(running)};
    }
}

long long findBestDelayPicoseconds(
    std::span<const long long> reference,
    std::span<const long long> target,
    long long coincWindowPs,
    long long delayStartPs,
    long long delayEndPs,
    long long delayStepPs,
    std::vector<std::pair<float, int>> *scratchResults) {
    std::vector<std::pair<float, int>> local;
    std::vector<std::pair<float, int>> &results =
        scratchResults ? *scratchResults : local;
    computeCoincidencesForRange(reference, target, coincWindowPs,
                                delayStartPs, delayEndPs, delayStepPs,
                                results);
    long long bestDelayPs = delayStartPs;
    int bestCount = std::numeric_limits<int>::min();
    for (const auto &entry : results) {
        const long long delayPs = static_cast<long long>(
            std::llround(static_cast<double>(entry.first) *
                         static_cast<double>(kPicosecondsPerNanosecond)));
        if (entry.second > bestCount) {
            bestCount = entry.second;
            bestDelayPs = delayPs;
        }
    }
    return bestDelayPs;
}

void writeResultsToFile(const std::vector<std::pair<float, int>> &results,
                        const std::string &filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    for (auto &p : results)
        out << p.first << "," << p.second << "\n";
    out.close();
}
