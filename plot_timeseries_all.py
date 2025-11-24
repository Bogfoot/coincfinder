"""
Time-series plots across all Delay_Scan_Data subfolders.

What it does
- Finds peak-delay for each "same" pair in the first folder only:
    HH (1,5), VV (2,6), DD (3,7), AA (4,8)   [DD/AA skipped in --setup4]
- For every folder (in sorted order) and for every second except the last one
  (to avoid partial seconds), it reads the coincidence counts at those fixed
  delays for the same pair and its corresponding cross pair:
    HV uses HH delay, VH uses VV delay; DA uses DD delay, AD uses AA delay.
- Builds continuous time axes by concatenating seconds from all folders.
- Plots:
    1) HV same+cross: HH, HV, VH, VV (shared X)
    2) DA same+cross: DD, DA, AD, AA (shared X; skipped for --setup4)
    3) Total visibility and total QBER (twin axes), where visibility/QBER are
       computed per-second from HV (and DA when available) same vs cross sums.

Usage (run next to CoincFinder outputs):
  python plot_timeseries_all.py
  python plot_timeseries_all.py --data-dir /path/to/Delay_Scan_Data
  python plot_timeseries_all.py --setup4
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Detector configuration
def detector_pairs(use_setup4: bool):
    if use_setup4:
        same_pairs = {"HH": (1, 4), "VV": (2, 3)}
        cross_pairs = {"HV": (1, 3), "VH": (2, 4)}
    else:
        same_pairs = {"HH": (1, 5), "VV": (2, 6), "DD": (3, 7), "AA": (4, 8)}
        cross_pairs = {"HV": (1, 6), "VH": (2, 5), "DA": (3, 8), "AD": (4, 7)}
    return same_pairs, cross_pairs

# CSV helpers
def parse_second(path: Path) -> int | None:
    name = path.stem  # delay_scan_<i>_vs_<j>_second_<sec>
    if "_second_" not in name:
        return None
    try:
        return int(name.split("_second_")[1])
    except ValueError:
        return None


def load_df(base: Path, i: int, j: int, sec: int):
    f = base / f"delay_scan_{i}_vs_{j}_second_{sec}.csv"
    if not f.exists():
        return None
    df = pd.read_csv(f, names=["delay_ns", "coinc"])
    return df if not df.empty else None

def first_second_peak_delay(folder: Path, pair: Tuple[int, int]) -> float | None:
    """
    Find the delay (ns) with max coincidences in second 0 only for a pair.
    """
    i, j = pair
    f = folder / f"delay_scan_{i}_vs_{j}_second_0.csv"
    if not f.exists():
        return None
    df = pd.read_csv(f, names=["delay_ns", "coinc"])
    if df.empty:
        return None
    idx = df["coinc"].idxmax()
    return float(df.loc[idx, "delay_ns"])

def count_at_delay(folder: Path, pair: Tuple[int, int], sec: int, delay: float) -> float:
    df = load_df(folder, pair[0], pair[1], sec)
    if df is None or df.empty:
        return 0.0
    idx = (df["delay_ns"] - delay).abs().idxmin()
    return float(df.loc[idx, "coinc"])


# Processing
def collect_delays(first_folder: Path, same_pairs: Dict[str, Tuple[int, int]]) -> Dict[str, float]:
    delays: Dict[str, float] = {}
    for label, pair in same_pairs.items():
        d = first_second_peak_delay(first_folder, pair)
        if d is not None:
            delays[label] = d
    return delays


def main():
    parser = argparse.ArgumentParser(description="Plot time-series coincidences across all delay-scan folders.")
    parser.add_argument("--data-dir", type=Path, default=Path("Delay_Scan_Data"))
    parser.add_argument("--setup4", action="store_true", help="Use 4-detector setup (HV only).")
    args = parser.parse_args()

    data_dir = args.data_dir
    folders = [p for p in sorted(data_dir.iterdir()) if p.is_dir()]
    if not folders:
        raise SystemExit(f"No folders found in {data_dir}")

    same_pairs, cross_pairs = detector_pairs(args.setup4)

    # Determine delays from first folder only
    delays = collect_delays(folders[0], same_pairs)
    if not delays:
        raise SystemExit("Could not determine peak delays in the first folder.")

    # Storage
    time_axis: List[int] = []
    counts: Dict[str, List[float]] = {k: [] for k in list(same_pairs.keys()) + list(cross_pairs.keys())}
    vis_total: List[float] = []
    qber_total: List[float] = []

    t = 0
    for folder in folders:
        seconds = sorted({parse_second(f) for f in folder.glob("delay_scan_*_second_*.csv") if parse_second(f) is not None})
        if len(seconds) <= 1:
            continue
        seconds = seconds[:-1]  # drop last second (potentially incomplete)

        for sec in seconds:
            # same counts
            same_vals = {}
            for label, pair in same_pairs.items():
                delay = delays.get(label)
                same_vals[label] = count_at_delay(folder, pair, sec, delay) if delay is not None else 0.0
                counts[label].append(same_vals[label])

            # cross counts
            cross_map = {
                "HV": ("HH", cross_pairs.get("HV")),
                "VH": ("VV", cross_pairs.get("VH")),
            }
            if not args.setup4:
                cross_map.update({
                    "DA": ("DD", cross_pairs.get("DA")),
                    "AD": ("AA", cross_pairs.get("AD")),
                })

            for clabel, (same_label, pair) in cross_map.items():
                if pair is None:
                    counts[clabel].append(0.0)
                    continue
                delay = delays.get(same_label)
                val = count_at_delay(folder, pair, sec, delay) if delay is not None else 0.0
                counts[clabel].append(val)

            # visibility/QBER per second
            hv_same = same_vals.get("HH", 0.0) + same_vals.get("VV", 0.0)
            hv_opp = counts["HV"][-1] + counts["VH"][-1]
            vis_hv = (hv_same - hv_opp) / (hv_same + hv_opp) if (hv_same + hv_opp) > 0 else np.nan
            q_hv = hv_opp / (hv_same + hv_opp) if (hv_same + hv_opp) > 0 else np.nan

            if args.setup4:
                vis_total.append(vis_hv)
                qber_total.append(q_hv)
            else:
                da_same = same_vals.get("DD", 0.0) + same_vals.get("AA", 0.0)
                da_opp = counts["DA"][-1] + counts["AD"][-1]
                vis_da = (da_same - da_opp) / (da_same + da_opp) if (da_same + da_opp) > 0 else np.nan
                q_da = da_opp / (da_same + da_opp) if (da_same + da_opp) > 0 else np.nan

                vis_total.append(np.nanmean([vis_hv, vis_da]))
                qber_total.append(np.nanmean([q_hv, q_da]))

            time_axis.append(t)
            t += 1

    # ---------------------------
    # Plots
    # ---------------------------
    if not time_axis:
        raise SystemExit("No data to plot (time axis is empty).")

    # HV plot
    fig_hv, ax_hv = plt.subplots(2, 1, sharex=True, figsize=(10, 6))
    ax_hv[0].plot(time_axis, counts["HH"], label="HH (same)", color="tab:blue")
    ax_hv[0].plot(time_axis, counts["VV"], label="VV (same)", color="tab:green")
    ax_hv[0].set_ylabel("Coincidences")
    ax_hv[0].set_title("HV basis: same")
    ax_hv[0].grid(True, alpha=0.3)
    ax_hv[0].legend()

    ax_hv[1].plot(time_axis, counts["HV"], label="HV (cross)", color="tab:orange")
    ax_hv[1].plot(time_axis, counts["VH"], label="VH (cross)", color="tab:red")
    ax_hv[1].set_ylabel("Coincidences")
    ax_hv[1].set_title("HV basis: cross")
    ax_hv[1].set_xlabel("Global seconds (concatenated)")
    ax_hv[1].grid(True, alpha=0.3)
    ax_hv[1].legend()

    # DA plot (only for 8-detector)
    if not args.setup4:
        fig_da, ax_da = plt.subplots(2, 1, sharex=True, figsize=(10, 6))
        ax_da[0].plot(time_axis, counts["DD"], label="DD (same)", color="tab:blue")
        ax_da[0].plot(time_axis, counts["AA"], label="AA (same)", color="tab:green")
        ax_da[0].set_ylabel("Coincidences")
        ax_da[0].set_title("DA basis: same")
        ax_da[0].grid(True, alpha=0.3)
        ax_da[0].legend()

        ax_da[1].plot(time_axis, counts["DA"], label="DA (cross)", color="tab:orange")
        ax_da[1].plot(time_axis, counts["AD"], label="AD (cross)", color="tab:red")
        ax_da[1].set_ylabel("Coincidences")
        ax_da[1].set_title("DA basis: cross")
        ax_da[1].set_xlabel("Global seconds (concatenated)")
        ax_da[1].grid(True, alpha=0.3)
        ax_da[1].legend()

    # Total visibility & QBER
    fig_tot, ax_tot = plt.subplots(figsize=(10, 4))
    ax_tot.plot(time_axis, np.array(vis_total) * 100, label="Total visibility", color="tab:blue")
    ax_tot.set_ylabel("Visibility (%)", color="tab:blue")
    ax_tot.tick_params(axis="y", labelcolor="tab:blue")
    ax_tot.set_xlabel("Global seconds (concatenated)")
    ax_tot.grid(True, alpha=0.3)

    ax_q = ax_tot.twinx()
    ax_q.plot(time_axis, np.array(qber_total) * 100, label="Total QBER", color="tab:red", linestyle="--")
    ax_q.set_ylabel("QBER (%)", color="tab:red")
    ax_q.tick_params(axis="y", labelcolor="tab:red")

    lines1, labels1 = ax_tot.get_legend_handles_labels()
    lines2, labels2 = ax_q.get_legend_handles_labels()
    ax_tot.legend(lines1 + lines2, labels1 + labels2, loc="upper right")
    ax_tot.set_title("Total visibility and QBER over time")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
