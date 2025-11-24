#pragma once

#include <map>
#include <vector>

#include "Singles.h"

/// Maintains per-channel Singles buckets for the last N seconds.
class RollingSingles {
public:
  explicit RollingSingles(long long windowSeconds = 200);

  /// Merge per-channel Singles produced by a chunk (e.g., readBINtoSingles).
  void appendChunk(const std::map<int, Singles> &chunk);

  /// Retrieve the Singles for `channel`. Returns empty instance when missing.
  const Singles &channelSingles(int channel) const;

  /// Latest chunk snapshots per channel (for histograms/auto-align).
  const std::vector<std::vector<Timestamp>> &latestChunk(int channel) const;

  /// Trim buckets older than the rolling window.
  void prune();

  /// Set the rolling window length in seconds.
  void setWindow(long long seconds);

  long long windowSeconds() const { return windowSeconds_; }
  long long latestSecond() const { return latestSecond_; }

  const std::map<int, Singles> &allChannels() const { return channels_; }

private:
  std::map<int, Singles> channels_;
  std::map<int, std::vector<std::vector<Timestamp>>> latestChunks_;
  long long windowSeconds_;
  long long latestSecond_;
};
