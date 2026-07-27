"""
Plot consecutive inter-arrival time histograms for every available channel.

Edit the values in the Config section below, then run:
  python3 plot_channel_consecutive_diffs.py
"""

from __future__ import annotations

from pathlib import Path

import coincfinder as cf
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Config
INPUT_FILE = Path(
    "../../CCode/TTDiffs/TimeStampsForCheckingCoincDistribution/Diffs.txt"
)
OUTPUT_DIR = Path("ChannelDiffs")
EXPOSURE_SECONDS = -1.0
MAX_BINS = 200
USE_LOG_X = True


def compute_consecutive_diffs(times_ps: np.ndarray) -> np.ndarray:
    if times_ps.size < 2:
        return np.empty(0, dtype=np.int64)
    return np.diff(times_ps)


def positive_bins(values: np.ndarray, max_bins: int) -> int:
    if values.size == 0:
        return 20
    return min(max_bins, max(20, values.size // 100))


def plot_channel_histogram(
    diffs_ps: np.ndarray,
    channel: int,
    out_path: Path,
    max_bins: int,
    log_x: bool,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))

    if diffs_ps.size == 0:
        ax.set_title(f"Channel {channel}: no consecutive differences available")
        ax.set_xlabel("Inter-arrival time (ps)")
        ax.set_ylabel("Count")
        ax.grid(alpha=0.25)
        fig.tight_layout()
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        return

    positive = diffs_ps[diffs_ps > 0]
    if log_x and positive.size:
        lo = positive.min()
        hi = positive.max()
        if lo == hi:
            bins = np.array([lo, hi + 1], dtype=float)
        else:
            bins = np.logspace(
                np.log10(lo), np.log10(hi), positive_bins(positive, max_bins)
            )
        ax.hist(positive, bins=bins, color="steelblue", edgecolor="black", alpha=0.7)
        ax.set_xscale("log")
        ax.set_xlabel("Inter-arrival time (ps, log scale)")
    else:
        ax.hist(
            diffs_ps,
            bins=positive_bins(diffs_ps, max_bins),
            color="steelblue",
            edgecolor="black",
            alpha=0.7,
        )
        ax.set_xlabel("Inter-arrival time (ps)")

    ax.set_title(f"Channel {channel}: consecutive arrival-time differences")
    ax.set_ylabel("Count")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> None:
    if not INPUT_FILE.exists():
        raise SystemExit(f"Input file not found: {INPUT_FILE}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    channels, duration_sec = cf.read_channels(str(INPUT_FILE), EXPOSURE_SECONDS)
    print(
        f"Loaded {INPUT_FILE} with duration {duration_sec:.3f} s "
        f"across channels {sorted(channels.keys())}"
    )

    summary_rows = []
    for channel in sorted(channels.keys()):
        times_ps = channels[channel]
        diffs_ps = compute_consecutive_diffs(times_ps)

        out_png = OUTPUT_DIR / f"channel_{channel}_consecutive_diffs.png"
        plot_channel_histogram(
            diffs_ps,
            channel=channel,
            out_path=out_png,
            max_bins=MAX_BINS,
            log_x=USE_LOG_X,
        )

        row = {
            "channel": channel,
            "events": int(times_ps.size),
            "diff_count": int(diffs_ps.size),
            "min_diff_ps": float(np.min(diffs_ps)) if diffs_ps.size else np.nan,
            "median_diff_ps": float(np.median(diffs_ps)) if diffs_ps.size else np.nan,
            "mean_diff_ps": float(np.mean(diffs_ps)) if diffs_ps.size else np.nan,
            "std_diff_ps": float(np.std(diffs_ps)) if diffs_ps.size else np.nan,
            "max_diff_ps": float(np.max(diffs_ps)) if diffs_ps.size else np.nan,
            "plot": out_png.name,
        }
        summary_rows.append(row)

        print(
            f"Channel {channel}: events={row['events']}, diffs={row['diff_count']}, "
            f"median={row['median_diff_ps']:.3f} ps, wrote {out_png}"
        )

    summary_path = OUTPUT_DIR / "channel_consecutive_diff_summary.csv"
    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    print(f"Wrote summary CSV to {summary_path}")


if __name__ == "__main__":
    main()
