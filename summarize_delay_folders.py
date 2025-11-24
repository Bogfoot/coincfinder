"""
Summarize CoincFinder delay-scan folders:
- Computes HV/DA visibility & QBER
- Computes total visibility = mean(HV_vis, DA_vis) when 8-detector
- Computes total QBER     = mean(HV_qber, DA_qber) when 8-detector
- Reports total coincidences per folder (HV, DA, overall)
- Plots:
    * HV & DA visibility together
    * HV & DA QBER together
    * Total visibility & total QBER on twin axes
    * Total coincidences per folder
- Saves summary table to summary.csv
Usage examples (run next to CoincFinder.exe or from project root):
  python summarize_delay_folders.py
  python summarize_delay_folders.py --data-dir C:/path/to/Delay_Scan_Data
  python summarize_delay_folders.py --setup4
"""

from __future__ import annotations

import argparse
from pathlib import Path
import os
import concurrent.futures
import itertools
import threading
import time
from typing import Iterable, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Detector setup
def detector_setup(use_setup4: bool):
    if use_setup4:
        pairs = {
            "H1–H2": (1, 4),
            "V1–V2": (2, 3),
        }
        same_HV = [(1, 4), (2, 3)]
        opp_HV = [(1, 3), (2, 4)]
        same_DA, opp_DA = [], []
    else:
        pairs = {
            "H1–H2": (1, 5),
            "V1–V2": (2, 6),
            "D1–D2": (3, 7),
            "A1–A2": (4, 8),
        }
        same_HV = [(1, 5), (2, 6)]
        opp_HV = [(1, 6), (2, 5)]
        same_DA = [(3, 7), (4, 8)]
        opp_DA = [(3, 8), (4, 7)]
    return pairs, same_HV, opp_HV, same_DA, opp_DA


# Helpers reading per-second CSVs
def parse_second_from_name(path: Path) -> int | None:
    name = path.stem  # delay_scan_<i>_vs_<j>_second_<sec>
    if "_second_" not in name:
        return None
    try:
        return int(name.split("_second_")[1])
    except ValueError:
        return None


def get_peak_delay_and_count(base: Path, i: int, j: int, sec: int):
    f = base / f"delay_scan_{i}_vs_{j}_second_{sec}.csv"
    if not f.exists():
        return None, 0.0
    df = pd.read_csv(f, names=["delay_ns", "coinc"])
    if df.empty:
        return None, 0.0
    idx_max = df["coinc"].idxmax()
    return float(df.loc[idx_max, "delay_ns"]), float(df.loc[idx_max, "coinc"])


def get_count_at_delay(base: Path, i: int, j: int, sec: int, delay_target: float):
    f = base / f"delay_scan_{i}_vs_{j}_second_{sec}.csv"
    if not f.exists():
        return 0.0
    df = pd.read_csv(f, names=["delay_ns", "coinc"])
    if df.empty:
        return 0.0
    idx = (df["delay_ns"] - delay_target).abs().idxmin()
    return float(df.loc[idx, "coinc"])


def counts_same_opp(base: Path, sec: int, same_pairs: Iterable[Tuple[int, int]], opp_pairs: Iterable[Tuple[int, int]]):
    C_same_total, C_opp_total = 0.0, 0.0
    for (si, sj), (oi, oj) in zip(same_pairs, opp_pairs):
        delay_peak, C_same = get_peak_delay_and_count(base, si, sj, sec)
        if delay_peak is None:
            continue
        C_opp = get_count_at_delay(base, oi, oj, sec, delay_peak)
        C_same_total += C_same
        C_opp_total += C_opp
    return C_same_total, C_opp_total


def compute_metrics(base: Path, seconds: Iterable[int], same_pairs, opp_pairs):
    vis, qber, same_counts, opp_counts = [], [], [], []
    for sec in seconds:
        C_same, C_opp = counts_same_opp(base, sec, same_pairs, opp_pairs)
        same_counts.append(C_same)
        opp_counts.append(C_opp)
        total = C_same + C_opp
        if total == 0:
            vis.append(np.nan)
            qber.append(np.nan)
        else:
            vis.append((C_same - C_opp) / total)
            qber.append(C_opp / total)
    same_arr = np.array(same_counts, dtype=float)
    opp_arr = np.array(opp_counts, dtype=float)
    return (
        np.array(vis, dtype=float),
        np.array(qber, dtype=float),
        same_arr,
        opp_arr,
        float(np.nansum(same_arr) + np.nansum(opp_arr)),
    )

# Main aggregation
def summarize_folder(folder: Path, use_setup4: bool):
    pairs, same_HV, opp_HV, same_DA, opp_DA = detector_setup(use_setup4)

    files = list(folder.glob("delay_scan_*_second_*.csv"))
    seconds = sorted({parse_second_from_name(f) for f in files if parse_second_from_name(f) is not None})
    if not seconds:
        return None  # no data

    vis_HV, q_HV, same_HV_counts, opp_HV_counts, total_HV = compute_metrics(folder, seconds, same_HV, opp_HV)

    vis_DA = q_DA = same_DA_counts = opp_DA_counts = np.array([])
    if not use_setup4:
        vis_DA, q_DA, same_DA_counts, opp_DA_counts, total_DA = compute_metrics(folder, seconds, same_DA, opp_DA)

    summary = {
        "folder": folder.name,
        "n_seconds": len(seconds),
        "HV_vis_mean": np.nanmean(vis_HV),
        "HV_vis_std": np.nanstd(vis_HV),
        "HV_qber_mean": np.nanmean(q_HV),
        "HV_qber_std": np.nanstd(q_HV),
        "HV_same_mean": np.nanmean(same_HV_counts),
        "HV_opp_mean": np.nanmean(opp_HV_counts),
        "HV_coinc_total": float(total_HV),
    }
    if not use_setup4:
        summary.update({
            "DA_vis_mean": np.nanmean(vis_DA),
            "DA_vis_std": np.nanstd(vis_DA),
            "DA_qber_mean": np.nanmean(q_DA),
            "DA_qber_std": np.nanstd(q_DA),
            "DA_same_mean": np.nanmean(same_DA_counts),
            "DA_opp_mean": np.nanmean(opp_DA_counts),
            "DA_coinc_total": float(total_DA),
        })
        # Total visibility/QBER as mean of HV & DA
        hv_vis = summary["HV_vis_mean"]
        da_vis = summary["DA_vis_mean"]
        hv_vis_std = summary["HV_vis_std"]
        da_vis_std = summary["DA_vis_std"]
        summary["Total_vis_mean"] = np.nanmean([hv_vis, da_vis])
        summary["Total_vis_std"] = np.sqrt((hv_vis_std**2 + da_vis_std**2) / 4.0)

        hv_q = summary["HV_qber_mean"]
        da_q = summary["DA_qber_mean"]
        hv_q_std = summary["HV_qber_std"]
        da_q_std = summary["DA_qber_std"]
        summary["Total_qber_mean"] = np.nanmean([hv_q, da_q])
        summary["Total_qber_std"] = np.sqrt((hv_q_std**2 + da_q_std**2) / 4.0)

        summary["Total_coinc_total"] = float(total_HV + total_DA)
    else:
        summary["Total_vis_mean"] = summary["HV_vis_mean"]
        summary["Total_vis_std"] = summary["HV_vis_std"]
        summary["Total_qber_mean"] = summary["HV_qber_mean"]
        summary["Total_qber_std"] = summary["HV_qber_std"]
        summary["Total_coinc_total"] = summary["HV_coinc_total"]
    return summary


def plot_series(ax, x: List[int], labels: List[str], values: List[float], title: str, ylabel: str, to_percent=False):
    vals = [v * 100 if to_percent else v for v in values]
    ax.plot(x, vals, marker="o", linestyle="-", color="tab:blue")
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.grid(True, alpha=0.3)


def main():
    parser = argparse.ArgumentParser(description="Summarize per-folder CoincFinder delay scans.")
    parser.add_argument("--data-dir", type=Path, default=Path("Delay_Scan_Data"), help="Path to Delay_Scan_Data root")
    parser.add_argument("--setup4", action="store_true", help="Use 4-detector setup instead of 8-detector")
    parser.add_argument("--no-plots", action="store_true", help="Skip plotting, just print table")
    parser.add_argument("--jobs", type=int, default=0, help="Parallel workers (0 = auto)")
    args = parser.parse_args()

    if not args.data_dir.exists():
        raise SystemExit(f"Data directory not found: {args.data_dir}")

    folders = [p for p in sorted(args.data_dir.iterdir()) if p.is_dir()]
    if not folders:
        raise SystemExit("No folders found in data directory.")

    # Parallelize per-folder summaries (I/O bound; threads are fine).
    workers = None if args.jobs <= 0 else args.jobs
    stop_spinner = threading.Event()

    def spinner():
        for ch in itertools.cycle("|/-\\"):
            if stop_spinner.is_set():
                break
            print(f"\rSummarizing folders... {ch}", end="", flush=True)
            time.sleep(0.1)
        # clear line
        print("\r" + " " * 40 + "\r", end="", flush=True)

    spinner_thread = threading.Thread(target=spinner, daemon=True)
    spinner_thread.start()

    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
            results = list(ex.map(lambda f: summarize_folder(f, args.setup4), folders))
    finally:
        stop_spinner.set()
        spinner_thread.join()

    summaries = [s for s in results if s]

    if not summaries:
        raise SystemExit("No delay scan data found.")

    df = pd.DataFrame(summaries)

    # Save summary CSV (comma-separated) into the data directory
    out_csv = args.data_dir / "summary.csv"
    df.to_csv(out_csv, index=False)
    print("\n=== Per-folder summary ===")
    print(df.to_string(index=False))
    print(f"\nSummary saved to: {out_csv}")

    if args.no_plots:
        return

    labels = df["folder"].tolist()
    x = list(range(len(labels)))

    if args.setup4:
        fig, axes = plt.subplots(1, 3, figsize=(16, 5))
        plot_series(axes[0], x, labels, df["HV_vis_mean"].tolist(), "HV Visibility", "Visibility (%)", to_percent=True)
        plot_series(axes[1], x, labels, df["HV_qber_mean"].tolist(), "HV QBER", "QBER (%)", to_percent=True)
        plot_series(axes[2], x, labels, df["Total_coinc_total"].tolist(), "Total coincidences", "Counts")
    else:
        fig, axes = plt.subplots(1, 4, figsize=(22, 5))
        plot_series(axes[0], x, labels, df["HV_vis_mean"].tolist(), "HV vs DA Visibility", "Visibility (%)", to_percent=True)
        axes[0].plot(x, np.array(df["DA_vis_mean"], float) * 100, marker="s", linestyle="--", color="tab:orange", label="DA visibility")
        axes[0].legend()

        plot_series(axes[1], x, labels, df["HV_qber_mean"].tolist(), "HV vs DA QBER", "QBER (%)", to_percent=True)
        axes[1].plot(x, np.array(df["DA_qber_mean"], float) * 100, marker="s", linestyle="--", color="tab:orange", label="DA QBER")
        axes[1].legend()

        # Total visibility & QBER on twin axes
        ax_tot = axes[2]
        ax_tot.plot(x, np.array(df["Total_vis_mean"], float) * 100, marker="o", linestyle="-", color="tab:blue", label="Total visibility")
        ax_tot.set_ylabel("Total visibility (%)", color="tab:blue")
        ax_tot.tick_params(axis="y", labelcolor="tab:blue")
        ax_tot.set_xticks(x)
        ax_tot.set_xticklabels(labels, rotation=45, ha="right")
        ax_tot.grid(True, alpha=0.3)
        ax_tot_twin = ax_tot.twinx()
        ax_tot_twin.plot(x, np.array(df["Total_qber_mean"], float) * 100, marker="s", linestyle="--", color="tab:red", label="Total QBER")
        ax_tot_twin.set_ylabel("Total QBER (%)", color="tab:red")
        ax_tot_twin.tick_params(axis="y", labelcolor="tab:red")
        ax_tot.set_title("Total visibility & QBER")
        # combine legends
        lines, labels_ = ax_tot.get_legend_handles_labels()
        lines2, labels2_ = ax_tot_twin.get_legend_handles_labels()
        ax_tot.legend(lines + lines2, labels_ + labels2_, loc="upper left")

        # Total coincidences
        plot_series(axes[3], x, labels, df["Total_coinc_total"].tolist(), "Total coincidences", "Counts")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
