# coincfinder

Fast coincidence scanning for time‑tagged detector singles and per‑pair timing analysis tools. The project ships:
- `CoincFinder` / `CoincPairs` CLIs for scanning and per-pair event dumping.
- `libcoincfinder` static library for reuse.
- Lightweight plotting/analysis scripts for delay sweeps and event timing inspection.

## Layout
- `src/` – CLI and core sources (`CoincFinder.cpp`, `Coincidences.cpp`, `ReadCSV.cpp`).
- `include/` – public headers mirroring the core sources.
- `build/` – out-of-source build (binaries, libraries, generated outputs).
- `Testing/` – quick sanity checks (`TestRolling.cpp`).
- Top-level helper scripts:
  - `plot_everything.py`, `plot_rates_all.py`, `plot_timeseries_all.py` – visualize scan outputs.
  - `Nicely_Plotted_various_diffs.py` – histogram/autocorr/phase-retrieval plots for event dumps.
  - `run_all_coincpairs.sh` – batch wrapper to run `CoincPairs` over a directory of `.bin` files.

## Prerequisites
- CMake ≥ 3.16
- C++20 compiler  
  - Linux: GCC ≥ 11 or Clang ≥ 13  
  - Windows: Visual Studio 2022 (MSVC)
- OpenMP (optional; auto-detected)
- Python ≥ 3.8 for helper scripts (matplotlib, pandas, numpy)

## Build
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```
Artifacts (Release):
- `build/CoincFinder`
- `build/CoincPairs`
- `build/libcoincfinder_core.a`

Clean rebuild if CMake cache tangles: `rm -rf build && cmake -S . -B build`.

## CoincFinder CLI
```
./CoincFinder <csv_or_bin> <coinc_window_ps> <delay_start_ns> <delay_end_ns> <delay_step_ns> <startSec> <stopSec>
```
Example:
```
./CoincFinder data/data.bin 250 8 12 0.01 0 600
```
Outputs go to `Delay_Scan_Data/` (per-second delay sweeps named `delay_scan_<ch1>_vs_<ch2>_second_<sec>.csv`).

## CoincPairs CLI (event dumps)
```
./CoincPairs <csv_or_bin> <coinc_window_ps> <delay_start_ns> <delay_end_ns> <delay_step_ns> <startSec> <stopSec> [rate_csv] --dump-events
```
- Emits per-pair event CSVs under `CoincEvents/<input_stem>/` when `--dump-events` is given.
- If `rate_csv` is provided, it is written inside that same per-run folder (see `run_all_coincpairs.sh`).
- See `docs/coinpairs.md` for a compact reference.

### Batch helper
`run_all_coincpairs.sh` runs `CoincPairs` over `8hMeasurement/*.bin`:
```bash
cd build
bash run_all_coincpairs.sh
```
Edit the variables at the top to change the coincidence window, delay range, or start/stop seconds. Each run writes its own `CoincEvents/<file_stem>/rate.csv` plus event dumps.

## Plotting event timing
After running `CoincPairs --dump-events`, use:
```bash
python3 Nicely_Plotted_various_diffs.py --event-dir CoincEvents --run <input_stem>
```
This tries to generate per-pair histograms, gap autocorr, and optional waveform reconstructions (phase retrieval and sqrt-FFT IFFT).

## Troubleshooting
- OpenMP missing: falls back to single-threaded; install an OpenMP-capable toolchain to speed up.
- Mixed configs: multi-config generators place outputs under `build/<Config>/`; use `--config Release` when running MSVC builds.
- Cache issues: delete `build/` and reconfigure.

## License
MIT — see `LICENSE`.
