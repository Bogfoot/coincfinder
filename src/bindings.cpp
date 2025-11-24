#include <cmath>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <span>

#include "Coincidences.h"
#include "ReadCSV.h"
#include "RollingSingles.h"
#include "Singles.h"

// Pybind11 module that mirrors the C++ CLI surface area. The bindings keep the
// docstrings short and defer to the underlying headers for deep detail, but the
// structure here explains how we translate Python types into spans and maps.

namespace py = pybind11;

PYBIND11_MODULE(coincfinder, m) {
  m.doc() = "Python bindings for the CoincFinder C++ library";

  // --- Bind Singles struct ---
  py::class_<Singles>(m, "Singles")
      .def(py::init<>())
      .def_readwrite("channel", &Singles::channel)
      .def_readwrite("base_second", &Singles::baseSecond)
      .def_readwrite("events_per_second", &Singles::eventsPerSecond)
      .def("__repr__", [](const Singles &s) {
        return "<Singles channel=" + std::to_string(s.channel) +
               ", seconds=" + std::to_string(s.eventsPerSecond.size()) + ">";
      });

  // --- Bind ReadCSV.h functions ---
  // Wrap duration out-parameters so Python gets a tuple (singles_map,
  // duration).
  m.def(
      "read_file_auto",
      [](const std::string &filename, double exposure_seconds) {
        double duration_sec = 0.0;
        auto singles = readFileAuto(filename, duration_sec, exposure_seconds);
        return std::make_pair(std::move(singles), duration_sec);
      },
      py::arg("filename"), py::arg("exposure_seconds") = -1.0,
      "Automatically read CSV or BIN file into a map<int, Singles>; returns "
      "(singles_map, measurement_duration_sec).");

  m.def(
      "read_csv_to_singles",
      [](const std::string &filename) {
        double duration_sec = 0.0;
        auto singles = readCSVtoSingles(filename, duration_sec);
        return std::make_pair(std::move(singles), duration_sec);
      },
      py::arg("filename"),
      "Read CSV file into map<int, Singles>; returns "
      "(singles_map, measurement_duration_sec).");

  m.def(
      "read_bin_to_singles",
      [](const std::string &filename) {
        double duration_sec = 0.0;
        auto singles = readBINtoSingles(filename, duration_sec);
        return std::make_pair(std::move(singles), duration_sec);
      },
      py::arg("filename"),
      "Read binary file into map<int, Singles>; returns "
      "(singles_map, measurement_duration_sec).");

  m.def("has_ending", &hasEnding, py::arg("string"), py::arg("ending"),
        "Check if a string ends with a given suffix");

  m.def("set_bucket_duration_seconds", &setBucketDurationSeconds,
        py::arg("seconds") = 1.0,
        "Set the time bucket duration (seconds) used when ingesting singles (default 1 s)."
        );
  m.def("get_bucket_duration_seconds", &bucketDurationSeconds,
        "Return the current bucket duration in seconds.");

  // --- Bind Coincidences.h functions ---
  // --- Count coincidences with delay (use ps everywhere in Python)
  m.def(
      "count_coincidences_with_delay_ps",
      [](const std::vector<long long> &ch1, const std::vector<long long> &ch2,
         double coinc_window_ps, double delay_ps) {
        const auto coinc_window_ll =
            static_cast<long long>(std::llround(coinc_window_ps));
        const auto delay_ll = static_cast<long long>(std::llround(delay_ps));
        return countCoincidencesWithDelay(std::span<const long long>(ch1),
                                          std::span<const long long>(ch2),
                                          coinc_window_ll, delay_ll);
      },
      py::arg("ch1"), py::arg("ch2"), py::arg("coinc_window_ps"),
      py::arg("delay_ps"), "Count coincidences (all arguments in picoseconds)");

  // --- Compute coincidences for range (accept delays in ps)
  m.def(
      "compute_coincidences_for_range_ps",
      [](const std::vector<long long> &ch1, const std::vector<long long> &ch2,
         double coinc_window_ps, double delay_start_ps, double delay_end_ps,
         double delay_step_ps) {
        std::vector<std::pair<float, int>> results;
        computeCoincidencesForRange(
            std::span<const long long>(ch1), std::span<const long long>(ch2),
            static_cast<long long>(std::llround(coinc_window_ps)),
            static_cast<long long>(std::llround(delay_start_ps)),
            static_cast<long long>(std::llround(delay_end_ps)),
            static_cast<long long>(std::llround(delay_step_ps)), results);
        return results;
      },
      py::arg("ch1"), py::arg("ch2"), py::arg("coinc_window_ps"),
      py::arg("delay_start_ps"), py::arg("delay_end_ps"),
      py::arg("delay_step_ps"),
      "Compute coincidences for delay range (all delays in picoseconds)");

  m.def(
      "compute_coincidences_for_range_hist_ps",
      [](const std::vector<long long> &ch1, const std::vector<long long> &ch2,
         double coinc_window_ps, double delay_start_ps, double delay_end_ps,
         double delay_step_ps) {
        std::vector<std::pair<float, int>> results;
        computeCoincidencesForRangeHistogram(
            std::span<const long long>(ch1), std::span<const long long>(ch2),
            static_cast<long long>(std::llround(coinc_window_ps)),
            static_cast<long long>(std::llround(delay_start_ps)),
            static_cast<long long>(std::llround(delay_end_ps)),
            static_cast<long long>(std::llround(delay_step_ps)), results);
        return results;
      },
      py::arg("ch1"), py::arg("ch2"), py::arg("coinc_window_ps"),
      py::arg("delay_start_ps"), py::arg("delay_end_ps"),
      py::arg("delay_step_ps"),
      "Histogram-based coincidence scan (all delays in picoseconds)");

  m.def(
      "count_nfold_coincidences",
      [](const std::vector<std::vector<long long>> &channels,
         double coinc_window_ps, const std::vector<long long> &offsets_ps) {
        std::vector<std::span<const long long>> spans;
        spans.reserve(channels.size());
        for (const auto &ch : channels)
          spans.emplace_back(ch.data(), ch.size());
        return countNFoldCoincidences(
            spans, static_cast<long long>(std::llround(coinc_window_ps)),
            std::span<const long long>(offsets_ps.data(), offsets_ps.size()));
      },
      py::arg("channels"), py::arg("coinc_window_ps"),
      py::arg("offsets_ps") = std::vector<long long>{},
      "Count N-fold coincidences across any number of channel traces "
      "(picoseconds)");

  m.def(
      "find_best_delay_ps",
      [](const std::vector<long long> &reference,
         const std::vector<long long> &target, double coinc_window_ps,
         double delay_start_ps, double delay_end_ps, double delay_step_ps) {
        return findBestDelayPicoseconds(
            std::span<const long long>(reference.data(), reference.size()),
            std::span<const long long>(target.data(), target.size()),
            static_cast<long long>(std::llround(coinc_window_ps)),
            static_cast<long long>(std::llround(delay_start_ps)),
            static_cast<long long>(std::llround(delay_end_ps)),
            static_cast<long long>(std::llround(delay_step_ps)));
      },
      py::arg("reference"), py::arg("target"), py::arg("coinc_window_ps"),
      py::arg("delay_start_ps"), py::arg("delay_end_ps"),
      py::arg("delay_step_ps"),
      "Return the delay (picoseconds) that maximizes coincidences between two "
      "channels.");

  py::class_<RollingSingles>(m, "RollingSingles")
      .def(py::init<long long>(), py::arg("window_seconds") = 200)
      .def("append_chunk", &RollingSingles::appendChunk, py::arg("chunk"))
      .def(
          "channel_singles",
          [](RollingSingles &self, int channel) -> const Singles & {
            return self.channelSingles(channel);
          },
          py::return_value_policy::reference_internal, py::arg("channel"))
      .def(
          "latest_chunk",
          [](RollingSingles &self,
             int channel) -> const std::vector<std::vector<Timestamp>> & {
            return self.latestChunk(channel);
          },
          py::return_value_policy::reference_internal, py::arg("channel"))
      .def("set_window_seconds", &RollingSingles::setWindow, py::arg("seconds"))
      .def("window_seconds", &RollingSingles::windowSeconds)
      .def("latest_second", &RollingSingles::latestSecond)
      .def(
          "all_channels",
          [](RollingSingles &self) -> const std::map<int, Singles> & {
            return self.allChannels();
          },
          py::return_value_policy::reference_internal);

  m.def("write_results_to_file", &writeResultsToFile, py::arg("results"),
        py::arg("filename"), "Write coincidence results to CSV");
}
