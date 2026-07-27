#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "Coincidences.h"

using Timestamp = long long;

namespace {
int naiveCoincidences(const std::vector<Timestamp> &ch1,
                      const std::vector<Timestamp> &ch2,
                      Timestamp windowPs,
                      Timestamp delayPs) {
    // windowPs is a full width (see Coincidences.h), so the acceptance range
    // is +/- half of it around the delay.
    const Timestamp halfWindow = windowPs / 2;
    const Timestamp lower = -halfWindow;
    const Timestamp upper = halfWindow;
    size_t i = 0;
    size_t j = 0;
    int count = 0;
    while (i < ch1.size() && j < ch2.size()) {
        const Timestamp shifted = ch1[i] - delayPs;
        const Timestamp diff = shifted - ch2[j];
        if (diff < lower) {
            ++i;
        } else if (diff > upper) {
            ++j;
        } else {
            ++count;
            ++i;
            ++j;
        }
    }
    return count;
}
} // namespace

void testHistogramMatchesNaive() {
    std::vector<Timestamp> ch1{0, 1'000, 2'000, 3'000, 4'000};
    std::vector<Timestamp> ch2{50, 1'050, 2'050, 3'050, 4'050};
    const Timestamp window = 100;
    const Timestamp delayStart = -200;
    const Timestamp delayEnd = 200;
    const Timestamp delayStep = 50;

    const DelayScan scan =
        scanDelays(ch1, ch2, window, delayStart, delayEnd, delayStep);
    assert(!scan.delaysPs.empty());
    for (size_t i = 0; i < scan.delaysPs.size(); ++i) {
        const Timestamp delayPs =
            static_cast<Timestamp>(std::llround(scan.delaysPs[i]));
        const int expected = naiveCoincidences(ch1, ch2, window, delayPs);
        assert(scan.counts[i] == expected);
    }
}

void testFindBestDelay() {
    std::vector<Timestamp> ref(30);
    for (size_t i = 0; i < ref.size(); ++i)
        ref[i] = static_cast<Timestamp>(i) * 2'000;
    const Timestamp offset = 1'250;
    std::vector<Timestamp> target(ref.size());
    for (size_t i = 0; i < ref.size(); ++i)
        target[i] = ref[i] + offset;

    // delay convention is ch1 - ch2 (see countCoincidencesWithDelay), so with
    // target = ref + offset every pair coincides for any delay in
    // [-offset-window/2, -offset+window/2] (window is a full width; constant
    // offset => constant plateau). findBestDelay isn't plateau-centered
    // (that's two-stage's job), so it returns the first grid point of that
    // tied-max plateau, i.e. its start.
    const Timestamp window = 200;
    const BestDelayResult result =
        findBestDelay(ref, target, window, -3'000, 3'000, 25);
    assert(result.bestDelayPs == -offset - window / 2);
}

void testFindBestDelayTwoStage() {
    // A clean, noise-free coincidence signal at a known offset: the coarse
    // scan should land near it and the fine scan should pin it down exactly.
    std::vector<Timestamp> ref(200);
    for (size_t i = 0; i < ref.size(); ++i)
        ref[i] = static_cast<Timestamp>(i) * 10'000;
    const Timestamp offset = 8'400;
    std::vector<Timestamp> target(ref.size());
    for (size_t i = 0; i < ref.size(); ++i)
        target[i] = ref[i] + offset;

    const TwoStageDelayResult result = findBestDelayTwoStage(
        ref, target, /*windowPs=*/200,
        /*coarseWindowPs=*/1'000, /*coarseHalfRangePs=*/50'000,
        /*coarseStepPs=*/500,
        /*fineHalfRangePs=*/2'000, /*fineStepPs=*/50);
    assert(result.bestDelayPs.has_value());
    assert(std::llround(*result.bestDelayPs) == -offset);
    assert(!result.fineScan.delaysPs.empty());
}

void testEstimateWindow() {
    // ch2 = ch1 - delay + jitter, jitter swept uniformly across [-halfJitter,
    // +halfJitter]. A uniform jitter distribution makes a flat-top histogram
    // of (ch1 - ch2) with width 2*halfJitter, so the FWHM (and thus the
    // estimated full window) should land close to 2*halfJitter.
    const int n = 400;
    const Timestamp delay = 5'000;
    const Timestamp halfJitter = 300;
    const Timestamp binPs = 20;
    std::vector<Timestamp> ch1(n), ch2(n);
    for (int i = 0; i < n; ++i) {
        ch1[i] = static_cast<Timestamp>(i) * 20'000;
        const Timestamp jitter =
            -halfJitter + (2 * halfJitter * i) / (n - 1);
        ch2[i] = ch1[i] - delay + jitter;
    }

    const WindowEstimate result = estimateCoincidenceWindow(
        ch1, ch2, delay, /*spanPs=*/2'000, binPs);
    assert(result.peakCount > result.backgroundCount);
    // A binPs-wide histogram can only resolve the true plateau edges to
    // within a couple of bins, so allow that much slack around 2*halfJitter.
    assert(std::abs(result.windowPs - static_cast<double>(2 * halfJitter)) <=
          3.0 * static_cast<double>(binPs));
}

void testNFoldCounts() {
    std::vector<Timestamp> base;
    for (size_t i = 0; i < 10; ++i)
        base.push_back(static_cast<Timestamp>(i) * 10'000);
    auto ch2 = base;
    auto ch3 = base;
    for (auto &ts : ch2)
        ts += 20;
    for (auto &ts : ch3)
        ts += 35;

    std::vector<std::span<const Timestamp>> spans = {
        std::span<const Timestamp>(base.data(), base.size()),
        std::span<const Timestamp>(ch2.data(), ch2.size()),
        std::span<const Timestamp>(ch3.data(), ch3.size()),
    };

    const int count = countNFoldCoincidences(spans, 100);
    assert(count == static_cast<int>(base.size()));

    spans.pop_back();
    const int pair = countNFoldCoincidences(spans, 100);
    assert(pair == static_cast<int>(base.size()));
}

int main() {
    testHistogramMatchesNaive();
    testFindBestDelay();
    testFindBestDelayTwoStage();
    testEstimateWindow();
    testNFoldCounts();
    std::cout << "All CoincFinder tests passed" << std::endl;
    return 0;
}
