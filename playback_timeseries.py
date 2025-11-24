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
    required = ["HH", "VV", "HV", "VH", "total_visibility", "total_qber"]
    for col in required:
        if col not in df.columns:
            raise SystemExit(f"timeseries.csv missing required column: {col}")
    has_da = {"DD", "AA", "DA", "AD"} <= set(df.columns)

    t_all = df["second"].to_numpy()
    data = {col: df[col].to_numpy() for col in df.columns if col != "second"}

    # Matplotlib setup (single figure with two rows)
    fig = plt.figure(figsize=(12, 8))
    gs = fig.add_gridspec(2, 1, height_ratios=[2, 1])
    ax_coinc = fig.add_subplot(gs[0, 0])

    # Combined coincidences (HV + DA if available)
    line_hh, = ax_coinc.plot([], [], label="HH", color="tab:blue")
    line_vv, = ax_coinc.plot([], [], label="VV", color="tab:green")
    line_hv, = ax_coinc.plot([], [], label="HV", color="tab:orange")
    line_vh, = ax_coinc.plot([], [], label="VH", color="tab:red")
    line_dd = line_aa = line_da = line_ad = None
    if has_da:
        line_dd, = ax_coinc.plot([], [], label="DD", color="tab:purple")
        line_aa, = ax_coinc.plot([], [], label="AA", color="tab:brown")
        line_da, = ax_coinc.plot([], [], label="DA", color="tab:pink")
        line_ad, = ax_coinc.plot([], [], label="AD", color="tab:gray")

    ax_coinc.set_ylabel("Coincidences")
    ax_coinc.set_title("Coincidences (same & cross)")
    ax_coinc.grid(True, alpha=0.3)
    ax_coinc.legend()

    # Total visibility/QBER (same figure, second row, twin y)
    ax_tot = fig.add_subplot(gs[1, 0])
    ax_q = ax_tot.twinx()
    line_vis, = ax_tot.plot([], [], label="Total visibility", color="tab:blue")
    line_q, = ax_q.plot([], [], label="Total QBER", color="tab:red", linestyle="--")
    ax_tot.set_ylabel("Visibility (%)", color="tab:blue")
    ax_q.set_ylabel("QBER (%)", color="tab:red")
    ax_tot.tick_params(axis="y", labelcolor="tab:blue")
    ax_q.tick_params(axis="y", labelcolor="tab:red")
    ax_tot.grid(True, alpha=0.3)
    ax_tot.set_xlabel("Global seconds")

    def update(frame):
        start = max(0, frame - args.window + 1)
        end = frame + 1
        if end <= start:
            return ()
        t = t_all[start:end]
        xlo, xhi = t[0] - 0.5, t[-1] + 0.5

        # Coincidences
        seg = {k: data[k][start:end] for k in data.keys()}
        line_hh.set_data(t, seg["HH"])
        line_vv.set_data(t, seg["VV"])
        line_hv.set_data(t, seg["HV"])
        line_vh.set_data(t, seg["VH"])
        if has_da:
            line_dd.set_data(t, seg["DD"])
            line_aa.set_data(t, seg["AA"])
            line_da.set_data(t, seg["DA"])
            line_ad.set_data(t, seg["AD"])

        y_candidates = [
            np.nanmax(seg["HH"]), np.nanmax(seg["VV"]),
            np.nanmax(seg["HV"]), np.nanmax(seg["VH"])
        ]
        if has_da:
            y_candidates += [
                np.nanmax(seg["DD"]), np.nanmax(seg["AA"]),
                np.nanmax(seg["DA"]), np.nanmax(seg["AD"])
            ]
        y_max = max([c for c in y_candidates if np.isfinite(c)] + [1])
        ax_coinc.set_xlim(xlo, xhi)
        ax_coinc.set_ylim(0, y_max * 1.1)
        ax_coinc.set_xlabel("Global seconds")

        # Update legend labels with latest values
        latest = {k: seg[k][-1] if len(seg[k]) else np.nan for k in seg}
        line_hh.set_label(f"HH ({latest.get('HH', float('nan')):.0f})")
        line_vv.set_label(f"VV ({latest.get('VV', float('nan')):.0f})")
        line_hv.set_label(f"HV ({latest.get('HV', float('nan')):.0f})")
        line_vh.set_label(f"VH ({latest.get('VH', float('nan')):.0f})")
        if has_da:
            line_dd.set_label(f"DD ({latest.get('DD', float('nan')):.0f})")
            line_aa.set_label(f"AA ({latest.get('AA', float('nan')):.0f})")
            line_da.set_label(f"DA ({latest.get('DA', float('nan')):.0f})")
            line_ad.set_label(f"AD ({latest.get('AD', float('nan')):.0f})")
        ax_coinc.legend(loc="upper left", ncol=2)

        # Total vis/QBER
        vis = seg["total_visibility"] * 100.0
        qber = seg["total_qber"] * 100.0
        line_vis.set_data(t, vis)
        line_q.set_data(t, qber)
        ax_tot.set_xlim(xlo, xhi)
        if len(vis):
            ax_tot.set_ylim(min(np.nanmin(vis), 0) - 5, max(np.nanmax(vis), 0) + 5)
            ax_q.set_ylim(min(np.nanmin(qber), 0) - 5, max(np.nanmax(qber), 0) + 5)

        # Update legend labels with latest vis/qber
        if len(vis):
            line_vis.set_label(f"Total visibility ({vis[-1]:.2f}%)")
        if len(qber):
            line_q.set_label(f"Total QBER ({qber[-1]:.2f}%)")
        lines1, labels1 = ax_tot.get_legend_handles_labels()
        lines2, labels2 = ax_q.get_legend_handles_labels()
        ax_tot.legend(lines1 + lines2, labels1 + labels2, loc="upper right")

        return ()

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
