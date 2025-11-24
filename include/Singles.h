#pragma once
#include <cstdint>
#include <cstddef>
#include <vector>

/// @file
/// Compact representation of time-tagged detector singles grouped into
/// contiguous one-second buckets. This structure is the backbone for both the
/// CLI and Python-facing APIs.

/// Alias for raw detector timestamps expressed in picoseconds.
using Timestamp = long long;

/// Represents singles collected on one detector channel, grouped by seconds.
struct Singles {
    /// Detector channel identifier (1-based).
    int channel = 0;
    /// Absolute second index associated with eventsPerSecond[0].
    long long baseSecond = 0;
    /// Per-second buckets of timestamps; bucket i => baseSecond + i.
    std::vector<std::vector<Timestamp>> eventsPerSecond;
};

/// Ensures the bucket for `second` exists and returns it for mutation.
inline std::vector<Timestamp> &ensureSecond(Singles &singles, long long second) {
    if (singles.eventsPerSecond.empty()) {
        singles.baseSecond = second;
        singles.eventsPerSecond.resize(1);
        return singles.eventsPerSecond.front();
    }

    if (second < singles.baseSecond) {
        const size_t prepend = static_cast<size_t>(singles.baseSecond - second);
        singles.eventsPerSecond.insert(singles.eventsPerSecond.begin(), prepend, {});
        singles.baseSecond = second;
        return singles.eventsPerSecond.front();
    }

    const long long maxSecond =
        singles.baseSecond + static_cast<long long>(singles.eventsPerSecond.size()) - 1;
    if (second > maxSecond) {
        const size_t newSize =
            static_cast<size_t>(second - singles.baseSecond + 1);
        singles.eventsPerSecond.resize(newSize);
    }

    const size_t idx = static_cast<size_t>(second - singles.baseSecond);
    return singles.eventsPerSecond[idx];
}

/// Returns the bucket for `second` or an empty view when out of range.
inline const std::vector<Timestamp> &eventsForSecond(const Singles &singles,
                                                     long long second) {
    static const std::vector<Timestamp> kEmpty;
    if (singles.eventsPerSecond.empty() || second < singles.baseSecond) {
        return kEmpty;
    }
    const size_t idx = static_cast<size_t>(second - singles.baseSecond);
    return idx < singles.eventsPerSecond.size() ? singles.eventsPerSecond[idx]
                                                : kEmpty;
}
