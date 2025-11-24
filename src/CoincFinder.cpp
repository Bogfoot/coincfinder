// CoincFinder CLI driver. Reads singles from CSV/BIN, scans a delay range for
// each detector pair, and writes per-second coincidence sweeps to disk.

#include <atomic>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <span>

#include "Coincidences.h"
#include "ReadCSV.h"
#include "Singles.h"

int main(int argc, char *argv[]) {
  if (argc < 8) {
    std::cerr
        << "Usage: " << argv[0]
        << " <csv_file> <coinc_window(ps)> <delay_start(ns)> <delay_end(ns)> "
           "<delay_step(ns)> <startSec> <stopSec>\n";
    return 1;
  }

  std::string csvFilename = argv[1];
  long long coincWindow = std::atoll(argv[2]);
  float delayStart = std::atof(argv[3]);
  float delayEnd = std::atof(argv[4]);
  float delayStep = std::atof(argv[5]);
  int startSec = std::atoi(argv[6]);
  int stopSec = std::atoi(argv[7]);

  if (coincWindow <= 0) {
    std::cerr << "Coincidence window must be positive.\n";
    return 1;
  }
  if (delayStep <= 0.0f) {
    std::cerr << "delay_step must be positive.\n";
    return 1;
  }
  if (delayEnd < delayStart) {
    std::cerr << "delay_end must be >= delay_start.\n";
    return 1;
  }
  if (startSec < 0 || stopSec < 0) {
    std::cerr << "startSec/stopSec must be non-negative.\n";
    return 1;
  }
  if (startSec > stopSec) {
    std::cerr << "startSec must be <= stopSec.\n";
    return 1;
  }

  // Convert once up front so the rest of the pipeline stays in integers.
  const auto nsToPs = [](float ns) -> long long {
    return static_cast<long long>(
        std::llround(static_cast<double>(ns) * 1000.0));
  };

  const long long delayStartPs = nsToPs(delayStart);
  const long long delayEndPs = nsToPs(delayEnd);
  const long long delayStepPs = nsToPs(delayStep);
  if (delayStepPs <= 0) {
    std::cerr << "delay_step too small once converted to picoseconds.\n";
    return 1;
  }

  std::cout << "Reading " << csvFilename << "...\n";
  double duration_sec = 0.0;
  auto singlesMap = readFileAuto(csvFilename, duration_sec);
  std::cout << "Measurement duration: " << duration_sec << " seconds\n";

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

  // Create folders
  try {
    std::filesystem::create_directories("Delay_Scan_Data");
  } catch (const std::exception &ex) {
    std::cerr << "Failed to create Delay_Scan_Data directory: " << ex.what()
              << "\n";
    return 1;
  }

  std::vector<std::pair<int, int>> coincidencePairs = {
      // Correlated (for Phi+)
      {1, 5}, // H-H
      {2, 6}, // V-V
      {3, 7}, // D-D
      {4, 8}, // A-A

      // Anti-correlated (for visibility)
      {1, 6}, // H-V
      {2, 5}, // V-H
      {3, 8}, // D-A
      {4, 7}  // A-D
  };

  // Build the subset of pairs that actually have data (avoids futile work).
  std::vector<std::pair<int, int>> activePairs;
  activePairs.reserve(coincidencePairs.size());
  for (const auto &pair : coincidencePairs) {
    if (singlesMap.count(pair.first) && singlesMap.count(pair.second)) {
      activePairs.push_back(pair);
    } else {
      std::cout << "Skipping ch" << pair.first << "-ch" << pair.second
                << " (missing singles).\n";
    }
  }
  if (activePairs.empty()) {
    std::cerr << "No coincidence pairs have data in the provided file.\n";
    return 1;
  }

  // Progress bar cuzz why not
  const int totalSeconds = stopSec - startSec + 1;
  const int totalJobs = static_cast<int>(activePairs.size()) * totalSeconds;
  std::atomic<int> jobsDone{0};

#pragma omp parallel for
  for (size_t p = 0; p < activePairs.size(); ++p) {
    const int ch1 = activePairs[p].first;
    const int ch2 = activePairs[p].second;
    const Singles &singles1 = singlesMap.at(ch1);
    const Singles &singles2 = singlesMap.at(ch2);

    std::vector<long long> mergedEvents;
    std::vector<std::pair<float, int>> results;
    size_t filesWritten = 0;
    for (int sec = startSec; sec <= stopSec; ++sec) {

      const auto &events1 = eventsForSecond(singles1, sec);
      if (events1.empty())
        continue;

      const auto &currentSecond = eventsForSecond(singles2, sec);
      const auto &nextSecond = eventsForSecond(singles2, sec + 1);
      if (currentSecond.empty() && nextSecond.empty())
        continue;

      // Include the first event from the next second so cross-second
      // coincidences survive
      const std::span<const long long> channel2Span =
          appendNextFirstEvent(currentSecond, nextSecond, mergedEvents);
      if (channel2Span.empty())
        continue;

      const std::span<const long long> channel1Span(events1.data(),
                                                    events1.size());

      std::string outFile = "Delay_Scan_Data/delay_scan_" +
                            std::to_string(ch1) + "_vs_" + std::to_string(ch2) +
                            "_second_" + std::to_string(sec) + ".csv";

      results.clear();
      computeCoincidencesForRange(channel1Span, channel2Span, coincWindow,
                                  delayStartPs, delayEndPs, delayStepPs,
                                  results);
      writeResultsToFile(results, outFile);
      ++filesWritten;

      int done = ++jobsDone;
      if (done == totalJobs || done % 50 == 0) {
#pragma omp critical
        std::cout << "\rProcessing " << done << " / " << totalJobs
                  << std::flush;
      }
    }

#pragma omp critical
    std::cout << "Finished ch" << ch1 << " vs ch" << ch2 << " (" << filesWritten
              << " seconds)\n";
  }

  if (totalJobs > 0) {
    std::cout << "\rProcessing " << totalJobs << " / " << totalJobs
              << " (done)\n";
  }

  std::cout << "\nSingles per second:\n";
  std::cout << "Second";
  for (int ch = 1; ch <= 8; ++ch)
    std::cout << "\tch" << ch;
  std::cout << "\n";

  long long maxSec = 0;
  for (auto &[ch, s] : singlesMap)
    if (!s.eventsPerSecond.empty()) {
      long long lastSecond =
          s.baseSecond + static_cast<long long>(s.eventsPerSecond.size()) - 1;
      maxSec = std::max(maxSec, lastSecond);
    }

  for (long long sec = 0; sec <= maxSec; ++sec) {
    std::cout << sec;
    for (int ch = 1; ch <= 8; ++ch) {
      long long count = 0;
      if (singlesMap.count(ch))
        count = eventsForSecond(singlesMap.at(ch), sec).size();
      std::cout << "\t" << count;
    }
    std::cout << "\n";
  }
  std::cout << "All done.\n";
  return 0;
}
