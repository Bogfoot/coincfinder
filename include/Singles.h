#pragma once
#include <algorithm>
#include <cstdint>
#include <cstddef>
#include <span>
#include <vector>

/// @file
/// Compact representation of time-tagged detector singles for one channel.
/// This structure is the backbone for both the CLI and Python-facing APIs.

/// Alias for raw detector timestamps expressed in picoseconds.
using Timestamp = long long;

/// Singles collected on one detector channel: a single flat, chronologically
/// sorted buffer (picoseconds, relative to the file's first timestamp), plus
/// the bucket width used to derive per-second slices on demand. Storing one
/// flat buffer instead of per-second buckets means there is nothing to
/// flatten/concatenate later - `events` already is the answer.
struct Singles {
    /// Detector channel identifier (1-based).
    int channel = 0;
    /// Width of one "second" bucket in picoseconds (set once at ingestion).
    long long bucketWidthPs = 0;
    /// Flat, ascending timestamps.
    std::vector<Timestamp> events;
};

/// Inserts `ts` keeping `singles.events` sorted. Real time-tag streams arrive
/// (near-)monotonically, so this is O(1) amortized in practice; only local
/// jitter pays for a shifted insert.
inline void insertSorted(Singles &singles, Timestamp ts) {
    auto &events = singles.events;
    if (events.empty() || ts >= events.back()) {
        events.push_back(ts);
        return;
    }
    events.insert(std::upper_bound(events.begin(), events.end(), ts), ts);
}

/// Returns the bucket index (elapsed seconds since the file's first
/// timestamp) that `ts` falls into.
inline long long secondOf(const Singles &singles, Timestamp ts) {
    return singles.bucketWidthPs > 0 ? ts / singles.bucketWidthPs : 0;
}

/// Returns the sub-span of `singles.events` covering bucket `second`, found
/// by binary search over the flat buffer (no per-second storage needed).
inline std::span<const Timestamp> eventsForSecond(const Singles &singles,
                                                   long long second) {
    if (singles.events.empty() || singles.bucketWidthPs <= 0)
        return {};

    const Timestamp lo = second * singles.bucketWidthPs;
    const Timestamp hi = lo + singles.bucketWidthPs;
    const auto begin = singles.events.begin();
    const auto end = singles.events.end();
    const auto first = std::lower_bound(begin, end, lo);
    const auto last = std::lower_bound(first, end, hi);

    const size_t startIdx = static_cast<size_t>(first - begin);
    const size_t count = static_cast<size_t>(last - first);
    return std::span<const Timestamp>(singles.events.data() + startIdx, count);
}

/// Bucket index of the earliest event, or 0 when empty.
inline long long firstSecond(const Singles &singles) {
    return singles.events.empty() ? 0 : secondOf(singles, singles.events.front());
}

/// Bucket index of the latest event, or -1 when empty.
inline long long lastSecond(const Singles &singles) {
    return singles.events.empty() ? -1 : secondOf(singles, singles.events.back());
}
