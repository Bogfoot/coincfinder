"""
Quick smoke test of the pybind11 coincfinder module using the current C++ lib.

What it does:
  - Loads the BIN file 2025-11-19_13_08_00_MDP_UVTP_exp_time_s_600.bin.
  - Reads singles via read_file_auto.
  - For each same pair (HH, VV, DD, AA) finds best delay, counts coincidences,
    and optionally collects timetag pairs via coincidences_with_delay_ps.
  - Reuses the same delays for cross pairs (HV/VH/DA/AD) and reports counts.
  - Writes a summary CSV and plots bar charts for counts and delays.
"""

from __future__ import annotations

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import coincfinder as cf

DATA_FILE = Path("2025-11-19_13_08_00_MDP_UVTP_exp_time_s_600.bin")
COINC_WINDOW_PS = 250
DELAY_START_NS = 8
DELAY_END_NS = 12
DELAY_STEP_NS = 0.01

SAME = {
    "HH": (1, 5),
    "VV": (2, 6),
    "DD": (3, 7),
    "AA": (4, 8),
}
CROSS = {
    "HV": (1, 6, "HH"),
    "VH": (2, 5, "VV"),
    "DA": (3, 8, "DD"),
    "AD": (4, 7, "AA"),
}


def main():
    if not DATA_FILE.exists():
        raise SystemExit(f"Data file not found: {DATA_FILE}")

    singles_map, duration = cf.read_file_auto(str(DATA_FILE))
    print(f"Loaded duration: {duration:.2f}s; channels: {list(singles_map.keys())}")

    def flatten(ch):
        v = singles_map[ch]
        flat = []
        for bucket in v.events_per_second:
            flat.extend(bucket)
        return np.array(flat, dtype=np.int64)

    delay_start_ps = int(DELAY_START_NS * 1000)
    delay_end_ps = int(DELAY_END_NS * 1000)
    delay_step_ps = int(DELAY_STEP_NS * 1000)

    summary = []

    delays_ns = {}
    for lbl, (c1, c2) in SAME.items():
        ch1 = flatten(c1)
        ch2 = flatten(c2)
        best = cf.find_best_delay_ps(ch1.tolist(), ch2.tolist(), COINC_WINDOW_PS,
                                     delay_start_ps, delay_end_ps, delay_step_ps)
        delays_ns[lbl] = best / 1000.0
        count = cf.count_coincidences_with_delay_ps(ch1.tolist(), ch2.tolist(),
                                                    COINC_WINDOW_PS, best)
        summary.append((lbl, delays_ns[lbl], count))

        # Collect a small sample of hits
        hits = cf.coincidences_with_delay_ps(ch1.tolist(), ch2.tolist(),
                                             COINC_WINDOW_PS, best, collect=True)
        pd.DataFrame(hits, columns=["t1_ps", "t2_ps"]).head(1000).to_csv(
            f"pybind_events_{lbl}.csv", index=False)

    # Cross using same delays
    for lbl, (c1, c2, base) in CROSS.items():
        if base not in delays_ns:
            continue
        delay_ps = int(delays_ns[base] * 1000)
        ch1 = flatten(c1)
        ch2 = flatten(c2)
        count = cf.count_coincidences_with_delay_ps(ch1.tolist(), ch2.tolist(),
                                                    COINC_WINDOW_PS, delay_ps)
        summary.append((lbl, delays_ns[base], count))

    df = pd.DataFrame(summary, columns=["pair", "delay_ns", "coinc"])
    df.to_csv("pybind_coinc_summary.csv", index=False)
    print(df)

    plt.figure(figsize=(8, 4))
    plt.bar(df["pair"], df["coinc"])
    plt.ylabel("Coincidences")
    plt.title("Counts at peak delays (pybind)")
    plt.tight_layout()
    plt.savefig("pybind_counts.png", dpi=150)

    plt.figure(figsize=(8, 4))
    plt.bar(df["pair"], df["delay_ns"])
    plt.ylabel("Delay (ns)")
    plt.title("Delays used per pair (pybind)")
    plt.tight_layout()
    plt.savefig("pybind_delays.png", dpi=150)

    print("Wrote pybind_coinc_summary.csv and plots.")


if __name__ == "__main__":
    main()
