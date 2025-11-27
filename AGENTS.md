# Repository Guidelines

## Project Structure & Modules
- Core sources live in `src/` (`CoincFinder.cpp` CLI, `Coincidences.cpp`, `ReadCSV.cpp`, `RollingSingles.cpp`).
- Public headers are in `include/` and mirror the core sources for library reuse.
- Lightweight tests live alongside code: `src/CoincFinderTests.cpp` (unit-style asserts) and `Testing/TestRolling.cpp` (integration sanity check).
- Helper analysis/plot scripts reside at the repo root (`plot_everything.py`, `playback_timeseries.py`, `plot_rates_all.py`, etc.); binary/CSV samples are stored beside them for quick experiments.
- Default out-of-source builds should go to `build/` (not checked in).

## Build, Run, and Development Commands
- Configure: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release` (add `-DOpenMP_CXX_FLAGS` if your toolchain needs help finding OpenMP).
- Build everything: `cmake --build build` (produces CLI `CoincFinder` and static lib `libcoincfinder.a`).
- Run CLI: `./build/CoincFinder <csv_or_bin> <coinc_window_ps> <delay_start_ns> <delay_end_ns> <delay_step_ns> <startSec> <stopSec>`.
- Clean rebuild when cache paths get tangled: `rm -rf build && cmake -S . -B build`.

## Testing Guidelines
- Unit assertions: after building, execute `./build/CoincFinderTests`; it must print “All CoincFinder tests passed”.
- Rolling buffer sanity check: `./build/Testing/TestRolling <sample.csv>` prints per-channel bucket counts; use it before large runs.
- Add new tests near the code they cover; keep them runnable via `ctest --test-dir build` if you register them in CMake.

## Coding Style & Naming
- C++20, 4-space indent, braces on the same line as declarations; prefer `std::` algorithms and `<span>` over raw pointers.
- File names are PascalCase with `.cpp/.h`; types use PascalCase, functions and variables use `camelCase`, constants may use `ALL_CAPS` if already established.
- Keep headers minimal and include what you use; avoid global state—pass containers by `const&` or `span`.
- Before pushing, format with `clang-format -style=Google` (matches existing layout) and build to catch warnings.

## Commit & Pull Request Guidelines
- Commit messages in history are short and emotive; move toward concise, imperative lines (e.g., `Fix delay clamping edge case`).
- Each PR should: describe the behavior change, list key commands run (`cmake --build build`, tests), and attach screenshots/plots if you changed output visualizations.
- Reference related issues or data files; note any new dependencies (OpenMP flags, Python packages).

## Security & Data Tips
- Sample time-tag data may be large; avoid committing new binaries—prefer links or small repro slices.
- Output directories like `Delay_Scan_Data/` and `build/` should stay untracked; add to `.gitignore` if new ones appear.
