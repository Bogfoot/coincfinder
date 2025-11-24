"""
Playback the aggregated time-series stored in Delay_Scan_Data/timeseries.csv.

Shows a rolling window (default 60 seconds) that slides forward every frame
(default frame interval 500 ms). Lines drop the oldest point and add the newest,
so you can watch how coincidences/visibility/QBER evolved across runs.

Run from the repo root or next to Delay_Scan_Data:
  python playback_timeseries.py
  python playback_timeseries.py --window 90 --interval 300 --data-dir ./Delay_Scan_Data
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Rolling playback of CoincFinder time-series.")
    parser.add_argument("--data-dir", type=Path, default=Path("Delay_Scan_Data"))
    parser.add_argument("--window", type=int, default=60, help="Window size in seconds (data samples) to display.")
    parser.add_argument("--interval", type=int, default=500, help="Frame interval in ms.")
    args = parser.parse_args()

    csv_path = args.data_dir / "timeseries.csv"
    if not csv_path.exists():
        raise SystemExit(f"timeseries.csv not found at {csv_path}. Run plot_timeseries_all.py first.")

    df = pd.read_csv(csv_path)
    if "second" not in df.columns:
        raise SystemExit("timeseries.csv missing 'second' column.")

    # Determine available channels
    has_da = {"DD", "AA", "DA", "AD"} <= set(df.columns)
    for col in ["HH", "VV", "HV", "VH", "total_visibility", "total_qber"]:
        if col not in df.columns:
            raise SystemExit(f"timeseries.csv missing required column: {col}")

    t_all = df["second"].to_numpy()
    data = {col: df[col].to_numpy() for col in df.columns if col != "second"}

    # Matplotlib setup
    fig = plt.figure(figsize=(12, 9))
    gs = fig.add_gridspec(3 if has_da else 2, 1, height_ratios=[1, 1, 1] if has_da else [1, 1])

    ax_hv = fig.add_subplot(gs[0, 0])
    ax_hv_cross = fig.add_subplot(gs[1, 0])
    ax_da = ax_da_cross = None
    if has_da:
        ax_da = fig.add_subplot(gs[2, 0])

    # HV same / cross
    line_hh, = ax_hv.plot([], [], label="HH (same)", color="tab:blue")
    line_vv, = ax_hv.plot([], [], label="VV (same)", color="tab:green")
    ax_hv.set_ylabel("Coincidences")
    ax_hv.set_title("HV basis (same)")
    ax_hv.legend()
    ax_hv.grid(True, alpha=0.3)

    line_hv, = ax_hv_cross.plot([], [], label="HV (cross)", color="tab:orange")
    line_vh, = ax_hv_cross.plot([], [], label="VH (cross)", color="tab:red")
    ax_hv_cross.set_ylabel("Coincidences")
    ax_hv_cross.set_title("HV basis (cross)")
    ax_hv_cross.grid(True, alpha=0.3)
    ax_hv_cross.legend()

    # DA same/cross if present
    if has_da:
        line_dd, = ax_da.plot([], [], label="DD (same)", color="tab:blue")
        line_aa, = ax_da.plot([], [], label="AA (same)", color="tab:green")
        line_da, = ax_da.plot([], [], label="DA (cross)", color="tab:orange")
        line_ad, = ax_da.plot([], [], label="AD (cross)", color="tab:red")
        ax_da.set_ylabel("Coincidences")
        ax_da.set_title("DA basis")
        ax_da.grid(True, alpha=0.3)
        ax_da.legend()

    # Total visibility/QBER twin axes
    fig_tot, ax_tot = plt.subplots(figsize=(10, 4))
    line_vis, = ax_tot.plot([], [], label="Total visibility", color="tab:blue")
    ax_tot.set_ylabel("Visibility (%)", color="tab:blue")
    ax_tot.tick_params(axis="y", labelcolor="tab:blue")
    ax_tot.grid(True, alpha=0.3)
    ax_tot.set_xlabel("Global seconds")

    ax_q = ax_tot.twinx()
    line_q, = ax_q.plot([], [], label="Total QBER", color="tab:red", linestyle="--")
    ax_q.set_ylabel("QBER (%)", color="tab:red")
    ax_q.tick_params(axis="y", labelcolor="tab:red")

    def update(frame):
        start = max(0, frame - args.window + 1)
        end = frame + 1
        t = t_all[start:end]

        # HV
        line_hh.set_data(t, data["HH"][start:end])
        line_vv.set_data(t, data["VV"][start:end])
        ax_hv.set_xlim(t[0], t[-1])
        y_max = max(np.nanmax(data["HH"][start:end]), np.nanmax(data["VV"][start:end]), 1)
        ax_hv.set_ylim(0, y_max * 1.1)

        line_hv.set_data(t, data["HV"][start:end])
        line_vh.set_data(t, data["VH"][start:end])
        ax_hv_cross.set_xlim(t[0], t[-1])
        y_max2 = max(np.nanmax(data["HV"][start:end]), np.nanmax(data["VH"][start:end]), 1)
        ax_hv_cross.set_ylim(0, y_max2 * 1.1)
        ax_hv_cross.set_xlabel("Global seconds")

        if has_da:
            line_dd.set_data(t, data["DD"][start:end])
            line_aa.set_data(t, data["AA"][start:end])
            line_da.set_data(t, data["DA"][start:end])
            line_ad.set_data(t, data["AD"][start:end])
            ax_da.set_xlim(t[0], t[-1])
            y_max3 = max(
                np.nanmax(data["DD"][start:end]),
                np.nanmax(data["AA"][start:end]),
                np.nanmax(data["DA"][start:end]),
                np.nanmax(data["AD"][start:end]),
                1,
            )
            ax_da.set_ylim(0, y_max3 * 1.1)
            ax_da.set_xlabel("Global seconds")

        # Total vis/QBER
        vis = data["total_visibility"][start:end] * 100.0
        qber = data["total_qber"][start:end] * 100.0
        line_vis.set_data(t, vis)
        line_q.set_data(t, qber)
        ax_tot.set_xlim(t[0], t[-1])
        if len(vis) > 0:
            ax_tot.set_ylim(min(np.nanmin(vis), 0) - 5, max(np.nanmax(vis), 0) + 5)
            ax_q.set_ylim(min(np.nanmin(qber), 0) - 5, max(np.nanmax(qber), 0) + 5)

        return (
            line_hh,
            line_vv,
            line_hv,
            line_vh,
            line_vis,
            line_q,
        )

    ani = animation.FuncAnimation(
        fig,
        update,
        frames=len(t_all),
        interval=args.interval,
        blit=False,
        repeat=False,
    )

    plt.show()


if __name__ == "__main__":
    main()
