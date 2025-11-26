"""
Compute best delays once, then apply them to all seconds to plot singles,
coincidences, visibility and QBER using the pybind `coincfinder` module.

Workflow:
  1) Load BIN/CSV via read_file_auto.
 2) Use a calibration second (default earliest) to find best delays for HH, VV, DD, AA (if present).
  3) For each second (optionally capped), count coincidences at fixed delays
     for same and cross pairs (HV uses HH delay, etc.).
  4) Plot:
       - Singles per second (sum across channels)
       - Coincidences per second (all pairs)
       - Visibility/QBER per second (HV and total if DA present)
  5) Save summary CSV: per-second counts, vis/QBER.

Usage:
  python plot_rates_all.py --file 2025-11-19_13_08_00_MDP_UVTP_exp_time_s_600.bin --seconds 600
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import coincfinder as cf

SAME = [("HH", 1, 5), ("VV", 2, 6), ("DD", 3, 7), ("AA", 4, 8)]
CROSS = [("HV", 1, 6, "HH"), ("VH", 2, 5, "VV"), ("DA", 3, 8, "DD"), ("AD", 4, 7, "AA")]
def main():
    ap = argparse.ArgumentParser(description="Fixed-delay per-second plots using coincfinder pybind module.")
    ap.add_argument("--file", type=Path, required=True, help="Input BIN/CSV file")
    ap.add_argument("--coinc-window-ps", type=float, default=250)
    ap.add_argument("--delay-start-ns", type=float, default=8)
    ap.add_argument("--delay-end-ns", type=float, default=12)
    ap.add_argument("--delay-step-ns", type=float, default=0.01)
    ap.add_argument("--seconds", type=int, default=None, help="Optional cap on number of seconds")
    ap.add_argument("--calib-second", type=int, default=None,
                    help="Absolute second used to calibrate delays (default: earliest available)")
    args = ap.parse_args()

    singles_map, duration = cf.read_file_auto(str(args.file))
    print(f"Loaded duration: {duration:.2f}s; channels: {list(singles_map.keys())}")

    delay_start_ps = int(args.delay_start_ns * 1000)
    delay_end_ps = int(args.delay_end_ns * 1000)
    delay_step_ps = int(args.delay_step_ns * 1000)

    # Determine available seconds
    available_seconds = set()
    for s in singles_map.values():
        available_seconds.update(range(s.base_second, s.base_second + len(s.events_per_second)))
    available_seconds = sorted(available_seconds)
    if not available_seconds:
        raise SystemExit("No data seconds found.")

    calib_sec = args.calib_second if args.calib_second is not None else available_seconds[0]
    if calib_sec not in available_seconds:
        raise SystemExit(f"Calibration second {calib_sec} not in data (available: {available_seconds[0]}..{available_seconds[-1]})")

    delays_ns = {}
    # Find best delays for same pairs
    for lbl, c1, c2 in SAME:
        if c1 not in singles_map or c2 not in singles_map:
            continue
        s1 = singles_map[c1]
        s2 = singles_map[c2]
        idx1 = calib_sec - s1.base_second
        idx2 = calib_sec - s2.base_second
        if not (0 <= idx1 < len(s1.events_per_second) and 0 <= idx2 < len(s2.events_per_second)):
            continue
        ch1 = s1.events_per_second[idx1]
        ch2 = s2.events_per_second[idx2]
        best = cf.find_best_delay_ps(ch1, ch2, args.coinc_window_ps,
                                     delay_start_ps, delay_end_ps, delay_step_ps)
        delays_ns[lbl] = best / 1000.0
        print(f"{lbl} best delay: {delays_ns[lbl]:.3f} ns (calib sec {calib_sec})")

    if not delays_ns:
        raise SystemExit("No delays found (missing channels?)")

    # Determine seconds range for processing
    seconds = available_seconds
    if args.seconds is not None:
        seconds = [sec for sec in seconds if sec - seconds[0] < args.seconds]

    rows = []
    singles_per_sec = {ch: [] for ch in singles_map.keys()}

    # Precompute which channels we actually need
    needed_channels = set(ch for _, ch1, ch2 in SAME for ch in (ch1, ch2))
    needed_channels.update(ch for _, ch1, ch2, _ in CROSS for ch in (ch1, ch2))

    for sec in seconds:
        row = {"second": sec}

        # Gather buckets once per second for needed channels
        buckets = {}
        for ch, s in singles_map.items():
            idx = sec - s.base_second
            if 0 <= idx < len(s.events_per_second):
                buckets[ch] = s.events_per_second[idx]
            else:
                buckets[ch] = []
            singles_per_sec[ch].append(len(buckets[ch]))

        # Same + cross pairs in one pass using a list
        all_pairs = []
        for lbl, c1, c2 in SAME:
            all_pairs.append((lbl, c1, c2, lbl))  # base delay = self
        for lbl, c1, c2, base in CROSS:
            all_pairs.append((lbl, c1, c2, base))

        for lbl, c1, c2, base in all_pairs:
            if base not in delays_ns or c1 not in buckets or c2 not in buckets:
                row[f"{lbl}_coinc"] = np.nan
                continue
            delay_ps = int(delays_ns[base] * 1000)
            ch1 = buckets[c1]
            ch2 = buckets[c2]
            if len(ch1) == 0 or len(ch2) == 0:
                row[f"{lbl}_coinc"] = np.nan
                continue
            row[f"{lbl}_coinc"] = cf.count_coincidences_with_delay_ps(
                ch1, ch2, args.coinc_window_ps, delay_ps)

        hh = row.get("HH_coinc") or 0
        vv = row.get("VV_coinc") or 0
        hv = row.get("HV_coinc") or 0
        vh = row.get("VH_coinc") or 0
        hv_same = hh + vv
        hv_opp = hv + vh
        hv_total = hv_same + hv_opp
        row["vis_hv"] = (hv_same - hv_opp) / hv_total if hv_total > 0 else np.nan
        row["qber_hv"] = hv_opp / hv_total if hv_total > 0 else np.nan

        dd = row.get("DD_coinc") or 0
        aa = row.get("AA_coinc") or 0
        da = row.get("DA_coinc") or 0
        ad = row.get("AD_coinc") or 0
        da_same = dd + aa
        da_opp = da + ad
        da_total = da_same + da_opp
        row["vis_da"] = (da_same - da_opp) / da_total if da_total > 0 else np.nan
        row["qber_da"] = da_opp / da_total if da_total > 0 else np.nan

        row["vis_total"] = np.nanmean([row["vis_hv"], row["vis_da"]])
        row["qber_total"] = np.nanmean([row["qber_hv"], row["qber_da"]])

        rows.append(row)

    df = pd.DataFrame(rows).sort_values("second").reset_index(drop=True)
    df.to_csv("per_second_fixed_delay_summary.csv", index=False)
    print("Saved per_second_fixed_delay_summary.csv")

    t = df["second"]

    # Singles per channel
    plt.figure(figsize=(10, 5))
    for ch, series in singles_per_sec.items():
        plt.plot(t, series, label=f"Ch{ch}")
    plt.xlabel("Second")
    plt.ylabel("Singles count")
    plt.title("Singles per second (per channel)")
    plt.grid(True, alpha=0.3)
    plt.legend(ncol=4, fontsize=8)
    plt.tight_layout()
    plt.savefig("singles_per_second.png", dpi=150)

    # Coincidences
    plt.figure(figsize=(10, 6))
    for col, label in [("HH_coinc", "HH"), ("VV_coinc", "VV"), ("HV_coinc", "HV"), ("VH_coinc", "VH"),
                       ("DD_coinc", "DD"), ("AA_coinc", "AA"), ("DA_coinc", "DA"), ("AD_coinc", "AD")]:
        if col in df:
            plt.plot(t, df[col], label=label)
    plt.xlabel("Second")
    plt.ylabel("Coincidences")
    plt.title("Per-second coincidences at fixed delays")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("coincidences_per_second.png", dpi=150)

    # Visibility/QBER
    plt.figure(figsize=(10, 4))
    plt.plot(t, df["vis_total"] * 100, label="Visibility", color="tab:blue")
    plt.plot(t, df["qber_total"] * 100, label="QBER", color="tab:red", linestyle="--")
    plt.xlabel("Second")
    plt.ylabel("Percent")
    plt.title("Visibility & QBER per second (fixed delays)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("vis_qber_per_second.png", dpi=150)

    print("Wrote singles_per_second.png, coincidences_per_second.png, vis_qber_per_second.png")


if __name__ == "__main__":
    main()
