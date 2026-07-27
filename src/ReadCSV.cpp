#include "ReadCSV.h"

// CSV/BIN ingestion helpers. They construct `Singles` objects in-place,
// inserting each timestamp directly into its channel's flat sorted buffer.

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <charconv>
#include <cmath>
#include <string_view>
#include <utility>

namespace {

constexpr long long kPicosecondsPerSecond = 1'000'000'000'000LL;
constexpr int kMaxChannels = 8;

std::atomic<double> gBucketSeconds{1.0};

template <typename T> bool parseIntegral(std::string_view token, T &value) {
  // Lightweight, locale-free parser so the loop stays allocation-free.
  const char *begin = token.data();
  const char *end = begin + token.size();
  while (begin < end && std::isspace(static_cast<unsigned char>(*begin))) {
    ++begin;
  }
  while (end > begin && std::isspace(static_cast<unsigned char>(*(end - 1)))) {
    --end;
  }
  if (begin == end) {
    return false;
  }

  auto result = std::from_chars(begin, end, value);
  return result.ec == std::errc() && result.ptr == end;
}

std::map<int, Singles>
finalizeSingles(std::array<Singles, kMaxChannels + 1> &channels) {
  std::map<int, Singles> result;
  for (int ch = 1; ch <= kMaxChannels; ++ch) {
    if (!channels[ch].events.empty()) {
      result.emplace(ch, std::move(channels[ch]));
    }
  }
  return result;
}

} // namespace

bool hasEnding(const std::string &str, const std::string &ending) {
  return str.size() >= ending.size() &&
         str.compare(str.size() - ending.size(), ending.size(), ending) == 0;
}

void setBucketDurationSeconds(double seconds) {
  const double clamped = seconds > 1e-9 ? seconds : 1.0;
  gBucketSeconds.store(clamped, std::memory_order_relaxed);
}

double bucketDurationSeconds() {
  return gBucketSeconds.load(std::memory_order_relaxed);
}

std::map<int, Singles> readFileAuto(const std::string &filename,
                                    double &duration_sec,
                                    double exposure_seconds) {
  if (exposure_seconds > 1e-9)
    setBucketDurationSeconds(exposure_seconds);
  if (hasEnding(filename, ".bin"))
    return readBINtoSingles(filename, duration_sec);
  return readCSVtoSingles(filename, duration_sec);
}

std::map<int, Singles> readCSVtoSingles(const std::string &filename,
                                        double &duration_sec) {
  std::ifstream file(filename);
  if (!file.is_open())
    throw std::runtime_error("Cannot open CSV file: " + filename);

  std::array<Singles, kMaxChannels + 1> channels;
  const long long bucketWidthPs = static_cast<long long>(
      std::llround(bucketDurationSeconds() * kPicosecondsPerSecond));
  for (int ch = 1; ch <= kMaxChannels; ++ch) {
    channels[ch].channel = ch;
    channels[ch].bucketWidthPs = bucketWidthPs;
  }

  std::string line;
  Timestamp firstTimestamp = 0;
  bool first = true;
  long long minTime = LLONG_MAX;
  long long maxTime = 0;

  while (std::getline(file, line)) {
    if (line.empty())
      continue;

    const size_t firstComma = line.find(',');
    if (firstComma == std::string::npos)
      continue;
    size_t secondComma = line.find(',', firstComma + 1);
    if (secondComma == std::string::npos)
      secondComma = line.size();

    std::string_view lineView(line);
    std::string_view timestampToken = lineView.substr(0, firstComma);
    std::string_view channelToken =
        lineView.substr(firstComma + 1, secondComma - firstComma - 1);

    Timestamp ts = 0;
    int ch = 0;
    if (!parseIntegral(timestampToken, ts) ||
        !parseIntegral(channelToken, ch)) {
      continue;
    }
    if (ch < 1 || ch > kMaxChannels || ts == 0)
      continue;

    if (first) {
      firstTimestamp = ts;
      first = false;
    }

    insertSorted(channels[ch], ts - firstTimestamp);

    if (ts < minTime)
      minTime = ts;
    if (ts > maxTime)
      maxTime = ts;
  }

  duration_sec = (maxTime > minTime) ? (maxTime - minTime) * 1e-12 : 0.0;
  return finalizeSingles(channels);
}

std::map<int, Singles> readBINtoSingles(const std::string &filename,
                                        double &duration_sec) {
  std::ifstream file(filename, std::ios::binary);
  if (!file.is_open())
    throw std::runtime_error("Cannot open BIN file: " + filename);

  file.seekg(40, std::ios::beg);

  std::array<Singles, kMaxChannels + 1> channels;
  const long long bucketWidthPs = static_cast<long long>(
      std::llround(bucketDurationSeconds() * kPicosecondsPerSecond));
  for (int ch = 1; ch <= kMaxChannels; ++ch) {
    channels[ch].channel = ch;
    channels[ch].bucketWidthPs = bucketWidthPs;
  }

  uint64_t t_raw = 0;
  uint16_t c_raw = 0;
  Timestamp firstTimestamp = 0;
  bool first = true;
  long long minTime = LLONG_MAX;
  long long maxTime = 0;

  while (file.read(reinterpret_cast<char *>(&t_raw), sizeof(t_raw))) {
    if (!file.read(reinterpret_cast<char *>(&c_raw), sizeof(c_raw)))
      break;

    Timestamp ts = static_cast<long long>(t_raw);
    int ch = static_cast<int>(c_raw) + 1;
    if (ch < 1 || ch > kMaxChannels || ts == 0)
      continue;

    if (first) {
      firstTimestamp = ts;
      first = false;
    }

    insertSorted(channels[ch], ts - firstTimestamp);

    if (ts < minTime)
      minTime = ts;
    if (ts > maxTime)
      maxTime = ts;
  }

  duration_sec = (maxTime > minTime) ? (maxTime - minTime) * 1e-12 : 0.0;
  return finalizeSingles(channels);
}
