#include "RollingSingles.h"

#include <algorithm>
#include <limits>

RollingSingles::RollingSingles(long long windowSeconds)
    : windowSeconds_(windowSeconds),
      latestSecond_(std::numeric_limits<long long>::min()) {}

void RollingSingles::appendChunk(const std::map<int, Singles> &chunk) {
    for (const auto &entry : chunk) {
        int channel = entry.first;
        const Singles &incoming = entry.second;
        if (incoming.eventsPerSecond.empty())
            continue;

        Singles &target = channels_[channel];
        if (target.eventsPerSecond.empty())
            target.channel = incoming.channel;

        auto &chunkSnapshot = latestChunks_[channel];
        chunkSnapshot = incoming.eventsPerSecond;

        for (size_t idx = 0; idx < incoming.eventsPerSecond.size(); ++idx) {
            long long second = incoming.baseSecond + static_cast<long long>(idx);
            latestSecond_ = std::max(latestSecond_, second);
            auto &bucket = ensureSecond(target, second);
            const auto &src = incoming.eventsPerSecond[idx];
            bucket.insert(bucket.end(), src.begin(), src.end());
        }
    }

    prune();
}

const Singles &RollingSingles::channelSingles(int channel) const {
    static const Singles kEmpty;
    auto it = channels_.find(channel);
    return it != channels_.end() ? it->second : kEmpty;
}

const std::vector<std::vector<Timestamp>> &RollingSingles::latestChunk(int channel) const {
    static const std::vector<std::vector<Timestamp>> kEmpty;
    auto it = latestChunks_.find(channel);
    return it != latestChunks_.end() ? it->second : kEmpty;
}

void RollingSingles::prune() {
    if (latestSecond_ == std::numeric_limits<long long>::min())
        return;
    const long long minSecond = latestSecond_ - windowSeconds_ + 1;
    for (auto &entry : channels_) {
        Singles &s = entry.second;
        if (s.eventsPerSecond.empty())
            continue;
        long long maxSecond = s.baseSecond + static_cast<long long>(s.eventsPerSecond.size()) - 1;
        if (maxSecond < minSecond)
            continue;
        if (s.baseSecond >= minSecond)
            continue;
        long long dropCount = minSecond - s.baseSecond;
        if (dropCount <= 0)
            continue;
        if (static_cast<size_t>(dropCount) >= s.eventsPerSecond.size()) {
            s.eventsPerSecond.clear();
            continue;
        }
        s.eventsPerSecond.erase(s.eventsPerSecond.begin(), s.eventsPerSecond.begin() + dropCount);
        s.baseSecond = minSecond;
    }
}

void RollingSingles::setWindow(long long seconds) {
    windowSeconds_ = std::max<long long>(1, seconds);
    prune();
}
