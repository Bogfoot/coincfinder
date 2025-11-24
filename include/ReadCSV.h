#pragma once
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <climits>

#include "Singles.h"

/// @file
/// Declares CSV/BIN ingestion helpers that populate `Singles` containers with
/// per-second buckets of timestamps. The accompanying implementation avoids
/// heap churn so coincidence scans can consume the data directly.

/// Configure bucket duration (seconds per time bucket). Defaults to 1s.
void setBucketDurationSeconds(double seconds);
double bucketDurationSeconds();

/// Dispatches to the appropriate reader based on filename suffix.
/// @param filename Path to CSV or BIN file.
/// @param duration_sec Populated with measurement duration in seconds.
std::map<int, Singles> readFileAuto(const std::string &filename, double &duration_sec,
                                    double exposure_seconds = -1.0);

/// Returns true if `str` ends with the requested suffix.
bool hasEnding(const std::string& str, const std::string& ending);

/// Parses a CSV file into per-channel singles (1-second buckets).
/// @param filename Path to CSV file (timestamp,channel,...).
/// @param duration Filled with measurement duration (seconds).
std::map<int, Singles> readCSVtoSingles(const std::string &filename, double &duration);

/// Parses a Qutools BIN file into per-channel singles.
std::map<int, Singles> readBINtoSingles(const std::string &filename, double &duration);
