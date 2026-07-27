#include <algorithm>
#include <cmath>
#include <optional>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <span>
#include <stdexcept>

#include "Coincidences.h"
#include "ReadCSV.h"
#include "Singles.h"

// Pybind11 module for the CoincFinder C++ library. All timing values (window,
// delay) are picoseconds throughout - there is only one unit, so none of the
// bindings carry a _ps/_np suffix. Array-like arguments (ch1, ch2, channels,
// offsets_ps) accept a Python list, tuple, or NumPy integer array via
// py::array::forcecast - one signature per operation.

namespace py = pybind11;

namespace {

using ArrayI64 = py::array_t<long long, py::array::c_style | py::array::forcecast>;

std::span<const long long> asSpan(const ArrayI64 &arr) {
  auto buf = arr.unchecked<1>();
  return std::span<const long long>(buf.data(0), static_cast<size_t>(buf.size()));
}

long long toPs(double valuePs) {
  return static_cast<long long>(std::llround(valuePs));
}

py::array_t<long long> toArray(std::span<const long long> v) {
  py::array_t<long long> out(static_cast<py::ssize_t>(v.size()));
  std::copy(v.begin(), v.end(), out.mutable_data());
  return out;
}

py::array_t<long long> toArray(const std::vector<long long> &v) {
  return toArray(std::span<const long long>(v));
}

py::array_t<double> toArray(const std::vector<double> &v) {
  py::array_t<double> out(static_cast<py::ssize_t>(v.size()));
  std::copy(v.begin(), v.end(), out.mutable_data());
  return out;
}

py::tuple scanToTuple(const DelayScan &scan) {
  return py::make_tuple(toArray(scan.delaysPs), toArray(scan.counts));
}

py::array_t<long long>
pairsToArray(const std::vector<std::pair<long long, long long>> &hits) {
  py::array_t<long long> out(
      {static_cast<py::ssize_t>(hits.size()), static_cast<py::ssize_t>(2)});
  auto buf = out.mutable_unchecked<2>();
  for (size_t i = 0; i < hits.size(); ++i) {
    buf(i, 0) = hits[i].first;
    buf(i, 1) = hits[i].second;
  }
  return out;
}

} // namespace

PYBIND11_MODULE(coincfinder, m) {
  m.doc() = "Python bindings for the CoincFinder C++ library";

  // --- Singles: per-channel flat, chronologically sorted timestamps (ps) ---
  py::class_<Singles>(m, "Singles")
      .def(py::init<>())
      .def_readonly("channel", &Singles::channel)
      .def_readonly("bucket_width_ps", &Singles::bucketWidthPs)
      .def_property_readonly(
          "events", [](const Singles &s) { return toArray(s.events); },
          "Flat, ascending int64 timestamp array (picoseconds).")
      .def(
          "events_for_second",
          [](const Singles &s, long long second) {
            return toArray(eventsForSecond(s, second));
          },
          py::arg("second"),
          "Timestamps (ps) falling within the given elapsed-second bucket.")
      .def("first_second", &firstSecond)
      .def("last_second", &lastSecond)
      .def("__repr__", [](const Singles &s) {
        return "<Singles channel=" + std::to_string(s.channel) +
               " events=" + std::to_string(s.events.size()) + ">";
      });

  // --- Reading ---
  m.def(
      "read_file_auto",
      [](const std::string &filename, double exposure_seconds) {
        double duration_sec = 0.0;
        auto singles = readFileAuto(filename, duration_sec, exposure_seconds);
        return std::make_pair(std::move(singles), duration_sec);
      },
      py::arg("filename"), py::arg("exposure_seconds") = -1.0,
      "Read CSV or BIN file into {channel: Singles}; returns "
      "(singles_map, measurement_duration_sec).");

  m.def(
      "read_csv_to_singles",
      [](const std::string &filename) {
        double duration_sec = 0.0;
        auto singles = readCSVtoSingles(filename, duration_sec);
        return std::make_pair(std::move(singles), duration_sec);
      },
      py::arg("filename"),
      "Read CSV file into {channel: Singles}; returns "
      "(singles_map, measurement_duration_sec).");

  m.def(
      "read_bin_to_singles",
      [](const std::string &filename) {
        double duration_sec = 0.0;
        auto singles = readBINtoSingles(filename, duration_sec);
        return std::make_pair(std::move(singles), duration_sec);
      },
      py::arg("filename"),
      "Read binary file into {channel: Singles}; returns "
      "(singles_map, measurement_duration_sec).");

  m.def(
      "read_channels",
      [](const std::string &filename, double exposure_seconds) {
        double duration_sec = 0.0;
        auto singlesMap = readFileAuto(filename, duration_sec, exposure_seconds);
        py::dict channels;
        for (auto &[ch, singles] : singlesMap)
          channels[py::int_(ch)] = toArray(singles.events);
        return py::make_tuple(channels, duration_sec);
      },
      py::arg("filename"), py::arg("exposure_seconds") = -1.0,
      "Read a CSV/BIN file straight into {channel: flat sorted int64 "
      "timestamp array (ps)}, plus measurement duration (s).");

  m.def(
      "read_channel",
      [](const std::string &filename, int channel, double exposure_seconds) {
        double duration_sec = 0.0;
        auto singlesMap = readFileAuto(filename, duration_sec, exposure_seconds);
        auto it = singlesMap.find(channel);
        if (it == singlesMap.end() || it->second.events.empty())
          throw std::runtime_error("channel " + std::to_string(channel) +
                                   " has no data in " + filename);
        return py::make_tuple(toArray(it->second.events), duration_sec);
      },
      py::arg("filename"), py::arg("channel"), py::arg("exposure_seconds") = -1.0,
      "Read a single channel as a flat int64 timestamp array (ps), plus "
      "measurement duration (s). Raises if the channel is absent or empty.");

  m.def("has_ending", &hasEnding, py::arg("string"), py::arg("ending"),
        "Check if a string ends with a given suffix");

  m.def("set_bucket_duration_seconds", &setBucketDurationSeconds,
        py::arg("seconds") = 1.0,
        "Set the time bucket duration (seconds) used when ingesting singles "
        "(default 1 s).");
  m.def("get_bucket_duration_seconds", &bucketDurationSeconds,
        "Return the current bucket duration in seconds.");

  // --- Pairwise coincidences ---
  m.def(
      "count_coincidences",
      [](ArrayI64 ch1, ArrayI64 ch2, double window_ps, double delay_ps) {
        return countCoincidencesWithDelay(asSpan(ch1), asSpan(ch2),
                                          toPs(window_ps), toPs(delay_ps));
      },
      py::arg("ch1"), py::arg("ch2"), py::arg("window_ps"), py::arg("delay_ps"),
      "Count coincidences between two channels at a fixed delay (picoseconds). "
      "window_ps is the full window width, centered on delay_ps.");

  m.def(
      "collect_coincidences",
      [](ArrayI64 ch1, ArrayI64 ch2, double window_ps, double delay_ps) {
        auto hits = collectCoincidencesWithDelay(asSpan(ch1), asSpan(ch2),
                                                  toPs(window_ps), toPs(delay_ps));
        return pairsToArray(hits);
      },
      py::arg("ch1"), py::arg("ch2"), py::arg("window_ps"), py::arg("delay_ps"),
      "Collect coincident timestamp pairs at a fixed delay; returns an "
      "(N, 2) int64 array of (t1_ps, t2_ps). window_ps is the full window "
      "width, centered on delay_ps.");

  m.def(
      "scan_delays",
      [](ArrayI64 ch1, ArrayI64 ch2, double window_ps, double delay_start_ps,
         double delay_end_ps, double delay_step_ps) {
        const auto scan =
            scanDelays(asSpan(ch1), asSpan(ch2), toPs(window_ps),
                      toPs(delay_start_ps), toPs(delay_end_ps),
                      toPs(delay_step_ps));
        return scanToTuple(scan);
      },
      py::arg("ch1"), py::arg("ch2"), py::arg("window_ps"),
      py::arg("delay_start_ps"), py::arg("delay_end_ps"),
      py::arg("delay_step_ps"),
      "Scan a delay range; returns (delays_ps, counts) as float64/int64 "
      "arrays. window_ps is the full window width, centered on each delay.");

  // --- Best delay ---
  m.def(
      "find_best_delay",
      [](ArrayI64 ch1, ArrayI64 ch2, double window_ps, double delay_start_ps,
         double delay_end_ps, double delay_step_ps) {
        const auto result =
            findBestDelay(asSpan(ch1), asSpan(ch2), toPs(window_ps),
                         toPs(delay_start_ps), toPs(delay_end_ps),
                         toPs(delay_step_ps));
        return py::make_tuple(static_cast<double>(result.bestDelayPs),
                              toArray(result.scan.delaysPs),
                              toArray(result.scan.counts));
      },
      py::arg("ch1"), py::arg("ch2"), py::arg("window_ps"),
      py::arg("delay_start_ps"), py::arg("delay_end_ps"),
      py::arg("delay_step_ps"),
      "Single-pass delay scan; returns (best_delay_ps, delays_ps, counts). "
      "window_ps is the full window width.");

  m.def(
      "find_best_delay_two_stage",
      [](ArrayI64 ch1, ArrayI64 ch2, double window_ps, double coarse_window_ps,
         double coarse_half_range_ps, std::optional<double> coarse_step_ps,
         double fine_half_range_ps, double fine_step_ps) {
        const double coarseStep = coarse_step_ps.value_or(coarse_window_ps / 2.0);
        const auto result = findBestDelayTwoStage(
            asSpan(ch1), asSpan(ch2), toPs(window_ps), toPs(coarse_window_ps),
            toPs(coarse_half_range_ps), toPs(coarseStep),
            toPs(fine_half_range_ps), toPs(fine_step_ps));
        py::object best =
            result.bestDelayPs ? py::cast(*result.bestDelayPs) : py::none();
        return py::make_tuple(best, toArray(result.fineScan.delaysPs),
                              toArray(result.fineScan.counts));
      },
      py::arg("ch1"), py::arg("ch2"), py::kw_only(), py::arg("window_ps"),
      py::arg("coarse_window_ps") = 1000.0,
      py::arg("coarse_half_range_ps") = 100'000.0,
      py::arg("coarse_step_ps") = py::none(),
      py::arg("fine_half_range_ps") = 50'000.0,
      py::arg("fine_step_ps") = 100.0,
      "Coarse scan over a wide range, then a fine scan + plateau-centered "
      "pick around the coarse peak. Returns (best_delay_ps or None, "
      "fine_delays_ps, fine_counts). window_ps/coarse_window_ps are full "
      "window widths.");

  m.def(
      "estimate_window",
      [](ArrayI64 ch1, ArrayI64 ch2, double delay_ps, double span_ps,
         double bin_ps) {
        const auto result = estimateCoincidenceWindow(
            asSpan(ch1), asSpan(ch2), toPs(delay_ps), toPs(span_ps),
            toPs(bin_ps));
        return py::make_tuple(result.windowPs, result.fwhmPs,
                              static_cast<double>(result.backgroundCount),
                              toArray(result.histogram.delaysPs),
                              toArray(result.histogram.counts));
      },
      py::arg("ch1"), py::arg("ch2"), py::arg("delay_ps"),
      py::arg("span_ps") = 5000.0, py::arg("bin_ps") = 20.0,
      "Estimate a coincidence window from the data: builds a fine histogram "
      "of (ch1 - ch2) around delay_ps and measures its FWHM against the "
      "surrounding accidental-coincidence background. Returns (window_ps, "
      "fwhm_ps, background_count, delays_ps, counts) - window_ps equals "
      "fwhm_ps (both full widths), a reasonable starting point for "
      "count_coincidences/scan_delays; compare background_count to "
      "counts.max() to judge how clean the peak is. Call "
      "find_best_delay_two_stage first to get delay_ps.");

  // --- N-fold ---
  m.def(
      "count_nfold_coincidences",
      [](py::list channels, double window_ps,
         std::optional<ArrayI64> offsets_ps) {
        std::vector<ArrayI64> owners;
        std::vector<std::span<const long long>> spans;
        owners.reserve(py::len(channels));
        spans.reserve(py::len(channels));
        for (auto item : channels) {
          owners.push_back(py::cast<ArrayI64>(item));
          spans.push_back(asSpan(owners.back()));
        }
        std::span<const long long> offsetsSpan;
        if (offsets_ps)
          offsetsSpan = asSpan(*offsets_ps);
        return countNFoldCoincidences(spans, toPs(window_ps), offsetsSpan);
      },
      py::arg("channels"), py::arg("window_ps"),
      py::arg("offsets_ps") = py::none(),
      "Count N-fold coincidences across any number of channels (picoseconds). "
      "channels is a list of array-like int64 timestamp arrays. window_ps is "
      "the full width within which all channels must have an event.");
}
