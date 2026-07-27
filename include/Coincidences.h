#pragma once
#include <fstream>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

/// @file
/// Core coincidence-accumulation helpers shared by both the CLI drivers and the
/// Python bindings. All timestamp/window/delay values are picoseconds,
/// passed around as lightweight spans to avoid redundant copies.
///
/// Delay convention: a hit is `ch1[i] - delayPs` landing within +/-
/// `coincWindowPs / 2` of `ch2[j]`, i.e. delayPs estimates `ch1 - ch2` at
/// coincidence. If ch2's events trail ch1's by `offset` (ch2 = ch1 + offset),
/// the best delay is `-offset`, not `+offset`.
///
/// Window convention: `coincWindowPs` (and `windowPs`/`binPs` below) is
/// always the *full* window width - a hit requires the difference to fall
/// within the window centered on the delay, i.e. +/- half of it.

/// Counts coincidences for a given delay between two channels (picoseconds).
int countCoincidencesWithDelay(std::span<const long long> ch1,
                               std::span<const long long> ch2,
                               long long coincWindowPs,
                               long long delayPs);

/// Collects timestamp pairs that fall within the coincidence window for a
/// given delay. Returns pairs of (t1_ps, t2_ps) in the original clock domain.
std::vector<std::pair<long long, long long>>
collectCoincidencesWithDelay(std::span<const long long> ch1,
                             std::span<const long long> ch2,
                             long long coincWindowPs,
                             long long delayPs);

/// Parallel-array result of a delay scan: delaysPs[i] pairs with counts[i].
struct DelayScan {
  std::vector<double> delaysPs;
  std::vector<long long> counts;
};

/// Scans a delay range and returns the (delay_ps, count) histogram using a
/// difference-array approach (single pass over the data).
DelayScan scanDelays(std::span<const long long> channel1,
                     std::span<const long long> channel2,
                     long long coincWindowPs,
                     long long delayStartPs,
                     long long delayEndPs,
                     long long delayStepPs);

/// Counts N-fold coincidences in a zero-delay window. When `channels.size()==2`
/// this simply calls `countCoincidencesWithDelay` with zero delay.
int countNFoldCoincidences(const std::vector<std::span<const long long>> &channels,
                           long long coincWindowPs,
                           std::span<const long long> offsetsPs = {});

/// Result of a single-pass best-delay search: the winning delay plus the full
/// scan histogram that produced it.
struct BestDelayResult {
  long long bestDelayPs = 0;
  DelayScan scan;
};

/// Finds the delay (picoseconds) within `[delayStartPs, delayEndPs]` that
/// yields the maximum coincidence count between `reference` and `target`.
BestDelayResult findBestDelay(std::span<const long long> reference,
                              std::span<const long long> target,
                              long long coincWindowPs,
                              long long delayStartPs,
                              long long delayEndPs,
                              long long delayStepPs);

/// Result of a coarse-then-fine best-delay search: `bestDelayPs` is unset when
/// neither stage found any coincidences.
struct TwoStageDelayResult {
  std::optional<double> bestDelayPs;
  DelayScan fineScan;
};

/// Coarse scan over `[-coarseHalfRangePs, +coarseHalfRangePs]`, plateau-centers
/// the coarse peak, then runs a fine scan +/- `fineHalfRangePs` around it at
/// `fineStepPs` resolution and plateau-centers again. Falls back to the coarse
/// pick if the fine scan finds nothing.
TwoStageDelayResult findBestDelayTwoStage(std::span<const long long> reference,
                                          std::span<const long long> target,
                                          long long windowPs,
                                          long long coarseWindowPs,
                                          long long coarseHalfRangePs,
                                          long long coarseStepPs,
                                          long long fineHalfRangePs,
                                          long long fineStepPs);

/// Estimate of the coincidence window suggested by the data's own timing
/// peak shape around a known delay.
struct WindowEstimate {
  double windowPs = 0.0;       // suggested full window width: == fwhmPs
  double fwhmPs = 0.0;         // full width at half maximum of the peak
  long long peakCount = 0;     // histogram count at the peak bin
  double backgroundCount = 0.0; // median per-bin count away from the peak
  DelayScan histogram;         // the fine histogram the estimate was read off
};

/// Builds a fine, near-unsmoothed histogram of `channel1 - channel2` around
/// `delayPs` (+/- `spanPs`, `binPs`-wide bins) and measures the FWHM of its
/// coincidence peak against the surrounding accidental-coincidence
/// background. `windowPs = fwhmPs` is a reasonable starting point for a
/// coincidence window (both are full widths); `backgroundCount` vs
/// `peakCount` indicates how noisy that estimate is.
WindowEstimate estimateCoincidenceWindow(std::span<const long long> channel1,
                                         std::span<const long long> channel2,
                                         long long delayPs,
                                         long long spanPs,
                                         long long binPs);

/// Writes a delay scan to `filename` as CSV (delay_ps,count).
void writeResultsToFile(const DelayScan &scan, const std::string &filename);

/// Returns a span over `currentSecond`, appending the first event from
/// `nextSecond` into `scratch` only when necessary. This mirrors the CLI logic
/// that preserves coincidences crossing a one-second boundary without copying
/// the entire bucket.
std::span<const long long>
appendNextFirstEvent(std::span<const long long> currentSecond,
                     std::span<const long long> nextSecond,
                     std::vector<long long> &scratch);
