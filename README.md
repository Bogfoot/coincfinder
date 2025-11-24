# coincfinder

A fast coincidence scanner for time‑tagged detector singles. The CLI ingests per‑channel timestamps (CSV or BIN), scans a user‑specified delay range, and writes per‑second coincidence sweeps to disk. A small static library (`coincfinder`) exposes the core logic for reuse in other applications.

## Layout
- `src/` – CLI (`CoincFinder.cpp`) and core sources (`Coincidences.cpp`, `ReadCSV.cpp`)
- `include/` – public headers (`Singles.h`, `Coincidences.h`, `ReadCSV.h`)
- `build/` – default out-of-source build output (executables, libraries, generated scan files)
- `plot_everything.py` – optional helper for visualizing results

## Prerequisites
- CMake ≥ 3.16
- C++20 compiler
  - Linux: GCC ≥ 11 or Clang ≥ 13
  - Windows: Visual Studio 2022 Build Tools (MSVC)
- OpenMP (optional) – if found, used for multithreading
- Python ≥ 3.8 (only for helper scripts; not required to build)

## Building

### Linux
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```
Artifacts land in `build/`:
- `CoincFinder` (CLI)
- `libcoincfinder.a` (static library; name may be `coincfinder.lib` on Windows)

### Windows (Visual Studio 2022 generator)
Run from a “x64 Native Tools Command Prompt for VS 2022”:
```cmd
cmake -S . -B build -G "Visual Studio 17 2022"
cmake --build build --config Release
```
Outputs are placed in `build/Release/` for the configuration build.

### Notes
- If you ever see cache/path mix-ups, delete the build directory and re-configure: `rm -rf build` (PowerShell/Bash) or `rmdir /s /q build` (cmd).
- On Windows the library target is named `coincfinder_lib` but its artifact is `coincfinder.lib`.

## CLI Usage
```
CoincFinder <csv_or_bin_file> <coinc_window_ps> <delay_start_ns> <delay_end_ns> <delay_step_ns> <startSec> <stopSec>
```

Example:
```
./CoincFinder data/data.bin 250 8 12 0.01 0 600
```
or
```
CoincFinder.exe data/data.bin 250 8 12 0.01 0 600
```

Parameters:
- `coinc_window_ps`  – coincidence window in picoseconds
- `delay_start_ns`   – starting delay (ns)
- `delay_end_ns`     – ending delay (ns)
- `delay_step_ns`    – step size (ns)
- `startSec/stopSec` – restrict processing to this second range

### Input formats
- `.bin` – binary time-tag format
- `.csv` – CSV time-tag format

## Output
- Results are written into `Delay_Scan_Data/` alongside the executable.
- For each detector pair and second, a file named `delay_scan_<ch1>_vs_<ch2>_second_<sec>.csv` is produced.

## Troubleshooting
- **Circular dependency / ResolveProjectReferences (MSBuild):** ensure the build directory is clean and reconfigure; name collisions are already avoided in `CMakeLists.txt`.
- **OpenMP not found:** build will fall back to single-threaded execution; install OpenMP-capable toolchain to enable.
- **Mixing debug/release artifacts:** multi-config generators place outputs under `build/<Config>/`; use `--config Release` when building and running.

## License
MIT — see `LICENSE` for details. Copyright (c) 2025 Adrian Udovičić, PhD student in physics at University of Ljubljana, Faculty of Mathematics and Physics.
