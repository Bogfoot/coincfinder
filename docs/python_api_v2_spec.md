# CoincFinder Python bindings: simplification spec

Status: proposal, not yet implemented. Written after prototyping a wrapper
against the current `.so` (see "Prototype" below) and auditing real call
sites across `Long-Distance-Entanglement-Distribution-FMF`, `QLaibLab`,
`LongDistanceQKD_1/2`, `TTSync`, `TTSyncBeyond`, and `Polarization
Controller`. Roughly 60 Python files import `coincfinder` across these
projects, so anything here should land additively, not as a breaking
rename (see "Rollout" at the end).

## 1. Confirmed bug: delay units flip between ps and ns

`computeCoincidencesForRange` (`src/Coincidences.cpp:175-233`) takes
`delayStartPs`/`delayEndPs`/`delayStepPs` in picoseconds, but writes the
delay column back out **divided by 1000**:

```cpp
// Coincidences.cpp:191-194
const long long delayPs =
    config.startPs + static_cast<long long>(idx) * config.stepPs;
results[idx] = {static_cast<float>(delayPs) / kPicosecondsPerNanosecond, 0};
```

So `compute_coincidences_for_range_ps`/`_np`/`_hist_*` take ps in and
return **nanoseconds** out, despite the `_ps` suffix and the docstring
("Compute coincidences for delay range (picoseconds)"). `find_best_delay_ps`
converts back to ps internally before returning
(`Coincidences.cpp:266-268`), so it's self-consistent, but the raw scan
functions are not.

This is exactly why `qkd_sync.py:343` has an unexplained
`* PS_PER_NS` after calling `compute_coincidences_for_range_np`. I hit the
same trap prototyping a fresh wrapper today: reasonable-looking code gave
a "best delay" of 10.6 ps instead of the correct 10,600 ps until I traced
it back to this line.

**Fix:** return the delay column in picoseconds, consistently, everywhere.
Drop the `/ kPicosecondsPerNanosecond` in `Coincidences.cpp:193` (and the
equivalent zero-channel early-return branch at the same spot).

## 2. Everything else wrong with the current shape

Evidence is `qkd_sync.py`, `test_pybind_api.py`, and
`QLaibLab/qlaiblib/io/coincfinder_backend.py` — three independent projects
that each rewrote the same workarounds:

- **Manual flattening, three times over.** `Singles.events_per_second` is
  `list[list[int]]` bucketed per second. Every caller writes the same
  `concatenate-then-sort` loop:
  `qkd_sync.py:flatten_channel`, `test_pybind_api.py:flatten`,
  `coincfinder_backend.py:_flatten_single`. This is ingestion-layer detail
  that shouldn't leak into Python at all.
- **Every op is tripled.** `count_coincidences_with_delay_{ps,np}`,
  `compute_coincidences_for_range_{ps,np,hist_ps,hist_np}`,
  `find_best_delay_{ps,np}` — six functions for three operations, split
  by whether the caller happens to have a `list` or an `ndarray`.
  pybind11's `py::array_t<int64_t, py::array::c_style |
  py::array::forcecast>` already accepts a Python list transparently (it
  forcecasts), so this split buys nothing. One signature per operation is
  sufficient.
- **No `_np` variant for collecting pairs.** `collect_coincidences_with_delay_ps`
  only takes `List[int]`, forcing `.tolist()` on every call
  (`qkd_sync.py:447-448`, `test_pybind_api.py:73-74`) — an O(n) Python-list
  copy of the largest arrays in the pipeline, only because the binding
  wasn't written to accept `py::array_t`.
- **The coarse→fine best-delay search is half-built in C++ and unused.**
  `findBestDelayPicoseconds` already accepts an optional
  `scratchResults` out-param for the histogram
  (`Coincidences.h:52-59`), but `bindings.cpp` never exposes it to Python.
  As a result, `qkd_sync.py` reimplements the entire two-stage
  coarse-then-fine search plus plateau-centering in ~90 lines
  (`find_best_delay`, `find_best_delay_near`,
  `best_delay_from_scan`, `qkd_sync.py:348-436`), and
  `QLaibLab/qlaiblib/coincidence/delays.py` does a version of the same
  thing again. Two separate projects independently rebuilt the same
  missing piece.
- **Scan results come back as `List[Tuple[float, int]]`.** Every caller
  immediately unzips it into two arrays
  (`qkd_sync.py:scan_delays:343-344`). Returning two NumPy arrays directly
  from the binding removes this loop everywhere it's called.
- **`std::pair<float, int>` loses precision.** The delay column is
  `float` (32-bit), i.e. ~7 significant digits, while picosecond
  timestamps run into the 10^12 range — a `float` can't represent a
  one-picosecond step at that magnitude without rounding. Should be
  `double`.

## 3. Prototype

`proto_wrapper.py` (pure-Python, wraps the existing `.so`, no C++ changes)
validates the shape below actually works end-to-end against
`coincfinder/2025-11-19_13_05_03_MDP_UVTP_exp_time_s_5.bin`, channels 1
and 5 (the "HH" pair from `test_pybind_api.py`):

```
duration=5.016s channels=[1..8]
best_delay_ps=10600.0   # matches test_pybind_api.py's known 8-12ns search band
count at best delay = 147
collected 147 timetag pairs, dtype=int64, shape=(147, 2)
```

That confirms both the target shape (flat channel dict, single
`scan_delays`/`find_best_delay` calls, ndarray in/out) and the unit bug
above (the prototype needed an explicit `* 1000.0` with a comment to get
a physically sane answer — that workaround should live in C++, not in
every Python project).

## 4. Proposed API

All functions below accept `array-like` (Python `list`, `tuple`, or any
NumPy integer array) via `py::array_t<int64_t, py::array::c_style |
py::array::forcecast>` — one binding, no `_ps`/`_np` split. All timing
values (windows, delays) are picoseconds, always; the name carries no
unit suffix because there is only one unit.

```python
# --- reading ---
def read_channels(filename: str, exposure_seconds: float = -1.0) \
        -> tuple[dict[int, np.ndarray], float]:
    """Flat, sorted int64 timestamp array per channel. Replaces
    read_file_auto + Singles + manual flatten. Internally still buckets
    by second for ingestion, but concatenates before returning."""

def read_channel(filename: str, channel: int, exposure_seconds: float = -1.0) \
        -> tuple[np.ndarray, float]:
    """Convenience: single channel, raises if absent/empty (mirrors the
    error-message quality of qkd_sync.py:load_channel)."""

# --- pairwise coincidences ---
def count_coincidences(ch1, ch2, window_ps: float, delay_ps: float) -> int: ...

def collect_coincidences(ch1, ch2, window_ps: float, delay_ps: float) \
        -> np.ndarray:  # shape (N, 2) int64, columns (t1_ps, t2_ps)
    ...

def scan_delays(ch1, ch2, window_ps: float,
                 delay_start_ps: float, delay_end_ps: float, delay_step_ps: float) \
        -> tuple[np.ndarray, np.ndarray]:  # (delays_ps float64, counts int64)
    """Fixes the ps-in/ns-out bug: delays_ps is picoseconds, matching the
    inputs, matching find_best_delay."""

# --- best delay ---
def find_best_delay(ch1, ch2, window_ps: float,
                     delay_start_ps: float, delay_end_ps: float, delay_step_ps: float) \
        -> tuple[float, np.ndarray, np.ndarray]:
    """(best_delay_ps, delays_ps, counts) - single-pass scan, histogram
    always returned (binds the scratchResults out-param that already
    exists in Coincidences.h but isn't exposed today)."""

def find_best_delay_two_stage(
    ch1, ch2, *, window_ps: float,
    coarse_window_ps: float = 1000.0,
    coarse_half_range_ps: float = 100_000.0,
    coarse_step_ps: float | None = None,   # default: coarse_window_ps / 2
    fine_half_range_ps: float = 50_000.0,
    fine_step_ps: float = 100.0,
) -> tuple[float | None, np.ndarray, np.ndarray]:
    """Coarse scan over a wide range, then a fine scan + plateau-centered
    pick around the coarse peak. Ports qkd_sync.py's find_best_delay +
    best_delay_from_scan into C++ so every project stops reimplementing
    it. Returns (best_delay_ps or None, fine_delays_ps, fine_counts)."""

# --- n-fold ---
def count_nfold_coincidences(channels: list, window_ps: float,
                              offsets_ps: list | None = None) -> int: ...
```

Internal C++ change needed to back `scan_delays`/`find_best_delay`
cleanly: change `std::vector<std::pair<float,int>>` to two parallel
`std::vector<double>`/`std::vector<long long>` (or a small struct), both
to fix the float-precision issue in §2 and because pybind11 can return
`py::array_t` for a `std::vector<double>` directly with a zero-copy cast,
avoiding the current `List[Tuple]` → Python-loop → two `np.ndarray`
unpacking dance every caller does today.

## 5. Naming

Drop `_ps` and `_np` suffixes on the new functions
(`count_coincidences`, `collect_coincidences`, `scan_delays`,
`find_best_delay`, `read_channels`, `count_nfold_coincidences`) — ps is
the only unit and array-like input is the only mode, so the suffix
encodes nothing. Keep `find_best_delay_two_stage`'s longer name since
it's doing something distinct from a single scan.

## 6. Rollout

~60 files across 6 sibling repos call the current API
(`Long-Distance-Entanglement-Distribution-FMF`, `QLaibLab`,
`LongDistanceQKD_1`, `LongDistanceQKD_2`, `TTSync`, `TTSyncBeyond`,
`Polarization Controller`). Two of them
(`qkd_sync.py`, `qlaiblib/coincidence/delays.py` +
`qlaiblib/io/coincfinder_backend.py`) already contain hand-rolled
versions of `find_best_delay_two_stage` and the flatten step — proof this
belongs in the library, not restated per project.

Recommend:
1. Add the new functions alongside the existing ones (this file's names
   are all new; nothing currently uses them, so no collision).
2. Fix the ps/ns bug in `computeCoincidencesForRange` itself (§1) — this
   changes the output of `compute_coincidences_for_range_{ps,np,hist_*}`
   for any *new* caller, but every *current* caller already multiplies by
   1000 to compensate (`qkd_sync.py:343`), or never uses the raw scan
   functions with fresh code. Grep each of the ~60 files for
   `compute_coincidences_for_range` before flipping this to be sure.
3. Migrate `qkd_sync.py` and `qlaiblib` to the new
   `find_best_delay_two_stage`/`read_channels` once built, deleting their
   local reimplementations.
4. Leave `read_file_auto`, `Singles`, and the `_ps`/`_np` pairs in place
   as deprecated but functional, rather than removing them immediately.

## 7. Aside: `.so` version drift (not in scope here, but worth knowing)

While prototyping, three different builds of `coincfinder.cpython-312-*.so`
were found active in this environment with different symbol sets:

- The repo-committed copy in `Long-Distance-Entanglement-Distribution-FMF/`
  has the `_np` variants and `collect_coincidences_with_delay_ps`, but no
  `RollingSingles`.
- `~/.local/lib/python3.12/site-packages/coincfinder...so` (what plain
  `pip`-based imports resolve to outside a project directory) has
  `RollingSingles` but *not* the `_np` variants, `collect_coincidences_with_delay_ps`,
  or `coincidences_with_delay_ps` — an older/different build entirely.
- The live source in `~/Documents/PhDCode/coincfinder/src/bindings.cpp`
  binds both `RollingSingles` and the `_np` variants, so it's newer than
  either installed copy.

Worth a single canonical build + install step (e.g. `pip install -e` from
the `coincfinder` source, removed from the other repos) once this spec is
implemented, so every project picks up the same binary.
