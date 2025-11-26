// CoincPairs CLI driver.
// Finds peak delays for same-channel pairs, then reports coincidence counts
// at those fixed delays for both same and cross pairs across the requested
// time window. Optionally dumps individual coincidence events (timetags).

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <span>
#include <string>
#include <optional>
#include <utility>
#include <vector>

#include "Coincidences.h"
#include "ReadCSV.h"
#include "Singles.h"

namespace {

struct PairInfo {
    int ch1;
    int ch2;
    std::string label;
    std::string delay_source; // which same-pair delay to reuse
};

struct DelayInfo {
    long long delayPs = 0;
    double delayNs = 0.0;
    bool valid = false;
};

std::span<const long long> spanWithNext(const Singles &s, int second,
                                        std::vector<long long> &scratch) {
    const auto &current = eventsForSecond(s, second);
    const auto &next = eventsForSecond(s, second + 1);
    return appendNextFirstEvent(current, next, scratch);
}

DelayInfo bestDelayForPair(const Singles &s1, const Singles &s2,
                           int second, long long coincWindowPs,
                           long long delayStartPs, long long delayEndPs,
                           long long delayStepPs,
                           std::vector<std::pair<float, int>> &scratchResults,
                           std::vector<long long> &scratch1,
                           std::vector<long long> &scratch2) {
    const auto span1 = spanWithNext(s1, second, scratch1);
    const auto span2 = spanWithNext(s2, second, scratch2);
    if (span1.empty() || span2.empty())
        return {};

    const long long delayPs =
        findBestDelayPicoseconds(span1, span2, coincWindowPs,
                                 delayStartPs, delayEndPs, delayStepPs,
                                 &scratchResults);
    DelayInfo info;
    info.delayPs = delayPs;
    info.delayNs = static_cast<double>(delayPs) / 1000.0;
    info.valid = true;
    return info;
}

int countAtDelay(const Singles &s1, const Singles &s2, int second,
                 long long coincWindowPs, long long delayPs,
                 std::vector<long long> &scratch1,
                 std::vector<long long> &scratch2) {
    const auto span1 = spanWithNext(s1, second, scratch1);
    const auto span2 = spanWithNext(s2, second, scratch2);
    if (span1.empty() || span2.empty())
        return 0;
    return countCoincidencesWithDelay(span1, span2, coincWindowPs, delayPs);
}

std::vector<std::pair<long long, long long>>
collectCoincidences(const Singles &s1, const Singles &s2, int second,
                    long long coincWindowPs, long long delayPs,
                    std::vector<long long> &scratch1,
                    std::vector<long long> &scratch2) {
    const auto span1 = spanWithNext(s1, second, scratch1);
    const auto span2 = spanWithNext(s2, second, scratch2);
    if (span1.empty() || span2.empty())
        return {};
    return collectCoincidencesWithDelay(span1, span2, coincWindowPs, delayPs);
}

} // namespace

void print_help(const char *exe) {
    std::cout
        << "CoincPairs - fixed-delay coincidence counter (optional timetags)\n"
        << "Usage: " << exe
        << " <csv|bin> <coinc_window_ps> <delay_start_ns> <delay_end_ns> <delay_step_ns> <startSec> <stopSec> [output_csv] [--dump-events]\n"
        << "Examples:\n"
        << "  " << exe << " data.bin 250 8 12 0.01 0 600\n"
        << "  " << exe << " data.bin 250 8 12 0.01 0 600 report.csv --dump-events\n\n"
        << "Behavior:\n"
        << "  - Finds peak delays for same pairs (HH, VV, DD, AA) at the first in-range second,\n"
        << "    reuses them for cross pairs (HV,VH,DA,AD).\n"
        << "  - Writes per-second counts to output_csv (default coincidences_report.csv).\n"
        << "  - With --dump-events, writes CoincEvents/<pair>.csv containing raw timetag pairs.\n"
        << "Notes:\n"
        << "  - startSec/stopSec are clamped to available data seconds.\n"
        << "  - delay_* in nanoseconds; window in picoseconds.\n";
}

int main(int argc, char *argv[]) {
    if (argc < 8) {
        print_help(argv[0]);
        return 1;
    }

    const std::string filename = argv[1];
    const long long coincWindowPs = std::atoll(argv[2]);
    const double delayStartNs = std::atof(argv[3]);
    const double delayEndNs = std::atof(argv[4]);
    const double delayStepNs = std::atof(argv[5]);
    int startSec = std::atoi(argv[6]);
    int stopSec = std::atoi(argv[7]);
    const bool dumpEvents =
        (argc >= 9 && std::string(argv[argc - 1]) == "--dump-events");

    const std::string outCsv =
        dumpEvents ? (argc >= 10 ? argv[8] : "coincidences_report.csv")
                   : (argc >= 9 ? argv[8] : "coincidences_report.csv");

    if (coincWindowPs <= 0 || delayStepNs <= 0.0 || delayEndNs < delayStartNs ||
        startSec < 0 || stopSec < 0 || startSec > stopSec) {
        std::cerr << "Invalid arguments.\n";
        return 1;
    }

    const auto nsToPs = [](double ns) {
        return static_cast<long long>(std::llround(ns * 1000.0));
    };
    const long long delayStartPs = nsToPs(delayStartNs);
    const long long delayEndPs = nsToPs(delayEndNs);
    const long long delayStepPs = nsToPs(delayStepNs);

    std::cout << "Reading " << filename << "...\n";
    double duration_sec = 0.0;
    auto singlesMap = readFileAuto(filename, duration_sec);

    long long earliestSec = std::numeric_limits<long long>::max();
    long long latestSec = std::numeric_limits<long long>::min();
    for (auto &[ch, singles] : singlesMap) {
        if (singles.eventsPerSecond.empty())
            continue;
        earliestSec = std::min(earliestSec, singles.baseSecond);
        long long last = singles.baseSecond +
                         static_cast<long long>(singles.eventsPerSecond.size()) - 1;
        latestSec = std::max(latestSec, last);
    }
    if (earliestSec == std::numeric_limits<long long>::max()) {
        std::cerr << "No singles data found.\n";
        return 1;
    }

    startSec = static_cast<int>(std::max<long long>(startSec, earliestSec));
    stopSec = static_cast<int>(std::min<long long>(stopSec, latestSec));
    if (startSec > stopSec) {
        std::cerr << "Requested second range has no overlap with data (available: "
                  << earliestSec << "-" << latestSec << ").\n";
        return 1;
    }

    // Define same and cross pairs
    std::vector<PairInfo> samePairs = {
        {1, 5, "HH", "HH"}, {2, 6, "VV", "VV"},
        {3, 7, "DD", "DD"}, {4, 8, "AA", "AA"}};
    std::vector<PairInfo> crossPairs = {
        {1, 6, "HV", "HH"}, {2, 5, "VH", "VV"},
        {3, 8, "DA", "DD"}, {4, 7, "AD", "AA"}};

    // Filter to existing channels
    auto hasChannel = [&](int ch) { return singlesMap.count(ch) > 0; };
    auto filterPairs = [&](std::vector<PairInfo> &pairs) {
        std::vector<PairInfo> filtered;
        for (auto &p : pairs) {
            if (hasChannel(p.ch1) && hasChannel(p.ch2))
                filtered.push_back(p);
        }
        pairs.swap(filtered);
    };
    filterPairs(samePairs);
    filterPairs(crossPairs);

    if (samePairs.empty()) {
        std::cerr << "No valid same-channel pairs found in data.\n";
        return 1;
    }

    // Compute best delays using the first available second in-range (startSec)
    std::map<std::string, DelayInfo> delays;
    std::vector<std::pair<float, int>> scratchResults;
    std::vector<long long> scratch1, scratch2;
    for (const auto &p : samePairs) {
        const Singles &s1 = singlesMap.at(p.ch1);
        const Singles &s2 = singlesMap.at(p.ch2);
        DelayInfo d = bestDelayForPair(s1, s2, startSec, coincWindowPs,
                                       delayStartPs, delayEndPs, delayStepPs,
                                       scratchResults, scratch1, scratch2);
        if (d.valid) {
            delays[p.label] = d;
            std::cout << "Delay " << p.label << ": " << d.delayNs << " ns\n";
        }
    }
    if (delays.empty()) {
        std::cerr << "Failed to determine any delays.\n";
        return 1;
    }

    // Prepare outputs
    std::ofstream out(outCsv);
    if (!out.is_open()) {
        std::cerr << "Cannot open output file: " << outCsv << "\n";
        return 1;
    }
    out << "second,pair,delay_ns,coincidences\n";

    std::map<std::string, std::ofstream> eventStreams;
    if (dumpEvents) {
        std::filesystem::create_directories("CoincEvents");
        auto openStream = [&](const std::string &label) {
            auto &stream = eventStreams[label];
            stream.open("CoincEvents/" + label + ".csv");
            stream << "second,t1_ps,t2_ps\n";
        };
        for (const auto &p : samePairs)
            openStream(p.label);
        for (const auto &p : crossPairs)
            openStream(p.label);
    }

    // Convenience list to process both same and cross with the same loop
    std::vector<PairInfo> allPairs;
    allPairs.reserve(samePairs.size() + crossPairs.size());
    allPairs.insert(allPairs.end(), samePairs.begin(), samePairs.end());
    allPairs.insert(allPairs.end(), crossPairs.begin(), crossPairs.end());

    // Process per-second counts
    for (int sec = startSec; sec <= stopSec; ++sec) {
        for (const auto &p : allPairs) {
            const auto itDelay = delays.find(p.delay_source);
            if (itDelay == delays.end() || !itDelay->second.valid)
                continue;

            const long long delayPs = itDelay->second.delayPs;
            const Singles &s1 = singlesMap.at(p.ch1);
            const Singles &s2 = singlesMap.at(p.ch2);

            const int count =
                countAtDelay(s1, s2, sec, coincWindowPs, delayPs, scratch1,
                             scratch2);
            out << sec << "," << p.label << "," << itDelay->second.delayNs
                << "," << count << "\n";

            if (dumpEvents) {
                auto hits = collectCoincidences(s1, s2, sec, coincWindowPs,
                                                delayPs, scratch1, scratch2);
                auto &stream = eventStreams[p.label];
                for (auto &h : hits) {
                    stream << sec << "," << h.first << "," << h.second << "\n";
                }
            }
        }
    }

    std::cout << "Wrote coincidence report to " << outCsv << "\n";
    if (dumpEvents) {
        std::cout << "Event CSVs written to CoincEvents/*.csv\n";
    }
    return 0;
}
