#pragma once
#include <fstream>
#include <span>
#include <string>
#include <utility>
#include <vector>

/// @file
/// Core coincidence-accumulation helpers shared by both the CLI driver and the
/// Python bindings. All timestamp/delay values are expressed in picoseconds and
/// passed around as lightweight spans to avoid redundant copies.

/// Counts coincidences (picoseconds) for a given delay between two channels.
int countCoincidencesWithDelay(std::span<const long long> ch1,
                               std::span<const long long> ch2,
                               long long coincWindowPs,
                               long long delayPs);

/// Scans a delay range and fills `results` with (delay_ns, coincidence_count)
/// using a histogram/difference-array approach (single pass over the data).
void computeCoincidencesForRange(std::span<const long long> channel1,
                                 std::span<const long long> channel2,
                                 long long coincWindowPs,
                                 long long delayStartPs, long long delayEndPs,
                                 long long delayStepPs,
                                 std::vector<std::pair<float, int>> &results);

inline void computeCoincidencesForRangeHistogram(
    std::span<const long long> channel1, std::span<const long long> channel2,
    long long coincWindowPs, long long delayStartPs, long long delayEndPs,
    long long delayStepPs, std::vector<std::pair<float, int>> &results) {
    computeCoincidencesForRange(channel1, channel2, coincWindowPs,
                                delayStartPs, delayEndPs, delayStepPs,
                                results);
}

/// Counts N-fold coincidences in a zero-delay window. When `channels.size()==2`
/// this simply calls `countCoincidencesWithDelay` with zero delay.
int countNFoldCoincidences(const std::vector<std::span<const long long>> &channels,
                           long long coincWindowPs,
                           std::span<const long long> offsetsPs = {});

/// Finds the delay (picoseconds) within `[delayStartPs, delayEndPs]` that yields
/// the maximum coincidence count between `reference` and `target`.
/// Returns the best delay in picoseconds and writes the histogram into
/// `scratchResults` when provided.
long long findBestDelayPicoseconds(
    std::span<const long long> reference,
    std::span<const long long> target,
    long long coincWindowPs,
    long long delayStartPs,
    long long delayEndPs,
    long long delayStepPs,
    std::vector<std::pair<float, int>> *scratchResults = nullptr);

/// Writes coincidence scan results to `filename` as CSV.
void writeResultsToFile(const std::vector<std::pair<float, int>> &results,
                        const std::string &filename);

/// Returns a span over `currentSecond`, appending the first event from
/// `nextSecond` into `scratch` only when necessary. This mirrors the CLI logic
/// that preserves coincidences crossing a one-second boundary without copying
/// the entire bucket.
std::span<const long long>
appendNextFirstEvent(const std::vector<long long> &currentSecond,
                     const std::vector<long long> &nextSecond,
                     std::vector<long long> &scratch);
