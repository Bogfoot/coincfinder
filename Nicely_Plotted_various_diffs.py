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
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

EVENT_DIR = Path("CoincEvents")
PAIRS = ["HH", "VV", "DD", "AA", "HV", "VH", "DA", "AD"]
# Channel mapping used by CoincPairs (first -> second)
PAIR_CHANNELS = {
    "HH": (1, 5),
    "VV": (2, 6),
    "DD": (3, 7),
    "AA": (4, 8),
    "HV": (1, 6),
    "VH": (2, 5),
    "DA": (3, 8),
    "AD": (4, 7),
}


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


def plot_autocorr(ax, data, title, max_lag):
    if data is None or len(data) < 2:
        ax.set_title(f"{title}\n(no data)")
        ax.axis("off")
        return
    x = data - np.mean(data)
    corr_full = np.correlate(x, x, mode="full")
    mid = len(corr_full) // 2
    corr = corr_full[mid : mid + max_lag + 1]
    corr = corr / corr[0] if corr[0] != 0 else corr
    lags = np.arange(len(corr))
    ax.plot(lags, corr, color="darkorange")
    ax.set_title(title)
    ax.set_xlabel("Lag (events)")
    ax.set_ylabel("Normalized autocorr")
    ax.grid(alpha=0.25)


def plot_rates(ax, seconds, counts, title):
    if len(seconds) == 0:
        ax.set_title(f"{title}\n(no data)")
        ax.axis("off")
        return
    ax.bar(seconds, counts, width=0.8, color="slategray")
    ax.set_title(title)
    ax.set_xlabel("Second")
    ax.set_ylabel("Coincidences")
    ax.grid(alpha=0.2)


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



def summarize_pair(pair, df, window_ps):
    t2t1 = df["t2_ps"].to_numpy() - df["t1_ps"].to_numpy()
    median = float(np.median(t2t1))
    mean = float(np.mean(t2t1))
    std = float(np.std(t2t1))
    mad = float(np.median(np.abs(t2t1 - median)))
    p16, p84 = np.percentile(t2t1, [16, 84])
    fwhm_like = p84 - p16
    summary = {
        "pair": pair,
        "median_ps": median,
        "mean_ps": mean,
        "std_ps": std,
        "mad_ps": mad,
        "p16_ps": float(p16),
        "p84_ps": float(p84),
        "fwhm_like_ps": float(fwhm_like),
    }
    if window_ps is not None:
        half = window_ps / 2.0
        frac_in_window = np.mean((t2t1 >= median - half) & (t2t1 <= median + half))
        summary["frac_in_window"] = float(frac_in_window)
    return summary


def process_pair(pair: str, window_ps: float | None, max_lag: int, channel_widths):
    path = EVENT_DIR / f"{pair}.csv"
    df = load_events(path)
    if df is None:
        print(f"{pair}: no events")
        return None

    diffs = compute_diffs(df)
    summary = summarize_pair(pair, df, window_ps)
    print(
        f"{pair}: median(t2-t1)={summary['median_ps']:.2f} ps, "
        f"FWHM-like={summary['fwhm_like_ps']:.2f} ps, "
        f"std={summary['std_ps']:.2f} ps"
        + (f", frac within window={summary['frac_in_window']*100:.1f}%"
           if window_ps is not None else "")
    )

    # Track channel jitter width for health hints
    if pair in PAIR_CHANNELS:
        ch1, ch2 = PAIR_CHANNELS[pair]
        for ch in (ch1, ch2):
            channel_widths.setdefault(ch, []).append(summary["fwhm_like_ps"])

    # Main histogram panel
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    fig.suptitle(f"{pair} coincidence timing distributions", fontsize=13, weight="bold")
    plot_hist(axes[0], diffs["distances"], f"{pair}: gap (next t1 - prev t2)", "gray")
    plot_hist(axes[1], diffs["coinc_times"], f"{pair}: t2 - t1", "steelblue")
    plot_hist(axes[2], diffs["consecutive_delta"], f"{pair}: Δ(t2-t1) between events", "crimson")
    plt.tight_layout()
    out_png = EVENT_DIR / f"{pair}_hist.png"
    plt.savefig(out_png, dpi=150)
    print(f"  wrote {out_png}")
    plt.show()
    plt.close(fig)

    # Autocorr + rate panel
    fig2, axes2 = plt.subplots(1, 2, figsize=(13, 4))
    plot_autocorr(axes2[0], diffs["distances"], f"{pair}: gap autocorr", max_lag)
    # Rates per second
    counts = df.groupby("second").size()
    plot_rates(axes2[1], counts.index.to_numpy(), counts.to_numpy(), f"{pair}: coincidences per second")
    plt.tight_layout()
    out_png2 = EVENT_DIR / f"{pair}_autocorr_rates.png"
    plt.savefig(out_png2, dpi=150)
    print(f"  wrote {out_png2}")
    plt.show()
    plt.close(fig2)

    return summary


def summarize_channels(channel_widths):
    if not channel_widths:
        print("Channel health: insufficient data.")
        return
    print("\nChannel health (average FWHM-like width across pairs, lower is better):")
    for ch, widths in sorted(channel_widths.items(), key=lambda kv: np.mean(kv[1])):
        print(f"  Ch {ch}: {np.mean(widths):.2f} ps over {len(widths)} pair(s)")


def main():
    parser = argparse.ArgumentParser(description="Plot CoincPairs event timing histograms.")
    parser.add_argument("--event-dir", type=Path, default=EVENT_DIR, help="Path to CoincEvents directory")
    parser.add_argument("--window-ps", type=float, default=None, help="Coincidence window used in CoincPairs (ps) for adequacy check")
    parser.add_argument("--max-lag", type=int, default=200, help="Max lag (events) for gap autocorrelation")
    args = parser.parse_args()

    event_dir = args.event_dir
    if not event_dir.exists():
        print(f"No CoincEvents directory found at {event_dir}")
        return
    available = [p for p in PAIRS if (event_dir / f"{p}.csv").exists()]
    if not available:
        print("No CoincEvents/<pair>.csv files found. Run CoincPairs with --dump-events first.")
        return

    global EVENT_DIR
    EVENT_DIR = event_dir  # so helper functions use the chosen dir

    summaries = []
    channel_widths = {}
    for pair in available:
        summary = process_pair(pair, args.window_ps, args.max_lag, channel_widths)
        if summary:
            summaries.append(summary)

    summarize_channels(channel_widths)


if __name__ == "__main__":
    main()
