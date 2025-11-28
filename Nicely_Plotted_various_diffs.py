"""
Plot coincidence timing deltas from CoincPairs event dumps.

CoincPairs writes raw events when invoked with --dump-events:
  CoincPairs <...> --dump-events
It produces CoincEvents/<PAIR>.csv with columns: second,t1_ps,t2_ps
for pairs: HH, VV, DD, AA, HV, VH, DA, AD (only those present in data).

This script reads the available CoincEvents/*.csv, computes a few
time-difference distributions, and saves per-pair PNGs.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

EVENT_DIR = Path("CoincEvents")
PAIRS = ["HH", "VV", "DD", "AA", "HV", "VH", "DA", "AD"]


def load_events(path: Path):
    if not path.exists():
        return None
    df = pd.read_csv(path)
    expected = {"second", "t1_ps", "t2_ps"}
    if not expected.issubset(df.columns):
        raise ValueError(f"{path} missing expected columns {expected}")
    if len(df) < 2:
        return None
    # Sort by (second, t1_ps) to ensure temporal order
    return df.sort_values(["second", "t1_ps", "t2_ps"], ignore_index=True)


def compute_diffs(df: pd.DataFrame):
    """Compute several timing deltas (all in picoseconds)."""
    t1 = df["t1_ps"].to_numpy()
    t2 = df["t2_ps"].to_numpy()

    diffs = {}
    # Gap between consecutive coincidences: next t1 minus previous t2
    diffs["distances"] = t1[1:] - t2[:-1]
    # Per-coincidence internal difference
    diffs["coinc_times"] = t2 - t1
    # Change in internal difference between consecutive coincidences
    diffs["consecutive_delta"] = (t2[1:] - t1[1:]) - (t2[:-1] - t1[:-1])
    return diffs


def plot_hist(ax, data, title, color):
    if data is None or len(data) == 0:
        ax.set_title(f"{title}\n(no data)")
        ax.axis("off")
        return
    n_bins = min(200, max(20, len(data) // 25))
    ax.hist(data, bins=n_bins, color=color, edgecolor="black", alpha=0.65)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("Time diff (ps)")
    ax.set_ylabel("Count")
    ax.grid(alpha=0.25)


def process_pair(pair: str):
    path = EVENT_DIR / f"{pair}.csv"
    df = load_events(path)
    if df is None:
        print(f"{pair}: no events")
        return

    diffs = compute_diffs(df)
    print(
        f"{pair}: mean(t2-t1)={np.nanmean(diffs['coinc_times']):.2f} ps, "
        f"median gap={np.nanmedian(diffs['distances']):.2f} ps"
    )

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    fig.suptitle(f"{pair} coincidence timing distributions", fontsize=13, weight="bold")
    plot_hist(axes[0], diffs["distances"], f"{pair}: gap (next t1 - prev t2)", "gray")
    plot_hist(axes[1], diffs["coinc_times"], f"{pair}: t2 - t1", "steelblue")
    plot_hist(axes[2], diffs["consecutive_delta"], f"{pair}: Δ(t2-t1) between events", "crimson")
    plt.tight_layout()
    out_png = EVENT_DIR / f"{pair}_hist.png"
    plt.savefig(out_png, dpi=150)
    print(f"  wrote {out_png}")
    # Show interactively for quick inspection; close after display so the next pair can render.
    plt.show()
    plt.close(fig)


def main():
    if not EVENT_DIR.exists():
        print(f"No CoincEvents directory found at {EVENT_DIR}")
        return
    available = [p for p in PAIRS if (EVENT_DIR / f"{p}.csv").exists()]
    if not available:
        print("No CoincEvents/<pair>.csv files found. Run CoincPairs with --dump-events first.")
        return
    for pair in available:
        process_pair(pair)


if __name__ == "__main__":
    main()
