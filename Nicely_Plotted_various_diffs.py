"""
Plot various time-difference distributions from pybind-dumped coincidence
events (files like pybind_events_<PAIR>.csv).
Each CSV is expected to have columns: t1_ps, t2_ps
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from pathlib import Path

PAIRS = ["HH", "VV", "DD", "AA", "HV", "VH", "DA", "AD"]
PREFIX = "pybind_events_"
SUFFIX = ".csv"

def load_events(pair):
    fname = Path(f"{PREFIX}{pair}{SUFFIX}")
    if not fname.exists():
        return None
    df = pd.read_csv(fname)
    if df.shape[0] < 2:
        return None
    return df

def compute_diffs(df):
    diffs = {}
    # distance between consecutive coincidences (t1 of i+1 minus t2 of i)
    diffs["distances"] = df["t1_ps"].iloc[1:].to_numpy() - df["t2_ps"].iloc[:-1].to_numpy()
    # coincidence time difference within each pair (t2 - t1)
    diffs["coinc_times"] = df["t2_ps"].to_numpy() - df["t1_ps"].to_numpy()
    # criss-cross: (next t2 - last t1) - (next t1 - last t2)
    diffs["criss_cross"] = (df["t2_ps"].iloc[1:].to_numpy() - df["t1_ps"].iloc[-1]).astype(float) - \
                           (df["t1_ps"].iloc[1:].to_numpy() - df["t2_ps"].iloc[-1]).astype(float)
    # consecutive coincidence time diffs
    diffs["consecutive"] = (df["t2_ps"].iloc[1:].to_numpy() - df["t1_ps"].iloc[1:].to_numpy()) - \
                           (df["t2_ps"].iloc[:-1].to_numpy() - df["t1_ps"].iloc[:-1].to_numpy())
    return diffs

def plot_hist(ax, data, title, color):
    if data is None or len(data)==0:
        ax.set_title(f"{title}\n(no data)")
        ax.axis("off")
        return
    n_bins = min(200, max(10, len(data)//20))
    ax.hist(data, bins=n_bins, color=color, edgecolor='black', alpha=0.6)
    ax.set_title(title, fontsize=12)
    ax.set_xlabel("Time diff (ps)")
    ax.set_ylabel("Count")
    ax.grid(alpha=0.3)

def main():
    available = [p for p in PAIRS if Path(f"{PREFIX}{p}{SUFFIX}").exists()]
    if not available:
        print("No pybind_events_*.csv files found.")
        return

    for pair in available:
        df = load_events(pair)
        if df is None:
            print(f"{pair}: no data")
            continue
        diffs = compute_diffs(df)
        print(f"{pair}: mean distances {np.nanmean(diffs['distances']) if len(diffs['distances']) else np.nan}")
        print(f"{pair}: mean coinc_times {np.nanmean(diffs['coinc_times']) if len(diffs['coinc_times']) else np.nan}")

        fig, axes = plt.subplots(2, 2, figsize=(12, 6))
        plot_hist(axes[0,0], diffs["distances"], f"{pair}: distances", 'gray')
        plot_hist(axes[0,1], diffs["coinc_times"], f"{pair}: t2-t1", 'blue')
        plot_hist(axes[1,0], diffs["criss_cross"], f"{pair}: criss-cross", 'green')
        plot_hist(axes[1,1], diffs["consecutive"], f"{pair}: consecutive Δ", 'red')
        plt.tight_layout()
        plt.savefig(f"pybind_events_{pair}_hist.png", dpi=150)
        plt.close(fig)

if __name__ == "__main__":
    main()
