# CoincPairs quick reference

CoincPairs scans time-tag files and optionally dumps raw per-pair coincidence events for offline timing analysis.

## Invocation
```
./CoincPairs <csv_or_bin> <coinc_window_ps> <delay_start_ns> <delay_end_ns> <delay_step_ns> <startSec> <stopSec> [rate_csv] --dump-events
```
- `coinc_window_ps` – coincidence half-window in picoseconds.
- `delay_start_ns` / `delay_end_ns` / `delay_step_ns` – delay sweep in nanoseconds.
- `startSec` / `stopSec` – restrict to a second range (0/0 processes full file).
- `rate_csv` (optional) – if provided, per-second singles/coincidence rates are written to this path.
- `--dump-events` – write per-pair event CSVs.

## Output layout
- `CoincEvents/<input_stem>/pair.csv` – columns: `second,t1_ps,t2_ps` for each available pair.
- `CoincEvents/<input_stem>/rate.csv` (if path supplied) – rate summary corresponding to that run.

## Batch runs
- Use `run_all_coincpairs.sh` (root) to process every `.bin` in `8hMeasurement/`; it places each run’s `rate.csv` inside its own `CoincEvents/<file_stem>/` folder.
- Adjust coincidence window, delay range, and start/stop seconds by editing the variables at the top of the script.

## Plotting dumps
After a run with `--dump-events`, inspect timing structure with:
```
python3 Nicely_Plotted_various_diffs.py --event-dir CoincEvents --run <input_stem>
```
This renders per-pair histograms, gap autocorrelation, and optional waveform reconstructions (phase retrieval and sqrt-FFT IFFT).
