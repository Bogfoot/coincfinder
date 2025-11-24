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
    const Timestamp lower = -windowPs;
    const Timestamp upper = windowPs;
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

    std::vector<std::pair<float, int>> results;
    computeCoincidencesForRange(ch1, ch2, window,
                                delayStart, delayEnd, delayStep,
                                results);
    assert(!results.empty());
    for (const auto &entry : results) {
        const Timestamp delayPs = static_cast<Timestamp>(
            std::llround(static_cast<double>(entry.first) * 1000.0));
        const int expected = naiveCoincidences(ch1, ch2, window, delayPs);
        assert(entry.second == expected);
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

    const Timestamp best = findBestDelayPicoseconds(ref, target,
                                                    200, -3'000, 3'000, 25);
    assert(best == offset);
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
    testNFoldCounts();
    std::cout << "All CoincFinder tests passed" << std::endl;
    return 0;
}
