#!/usr/bin/env python3
"""
Concatenate and plot every CoincEvents/*/rate.csv back-to-back.

Usage: adjust the variables in the “Config” section below, then run
  python plot_rate_consecutive.py
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def load_rate_tables(data_dir: Path) -> Tuple[pd.DataFrame, List[int]]:
    """
    Load every rate.csv under data_dir/*, adding an absolute second column
    that keeps time continuous across runs. Also returns boundaries (seconds)
    where runs change for optional vertical markers.
    """
    folders = sorted(p for p in data_dir.iterdir() if p.is_dir())
    frames: List[pd.DataFrame] = []
    boundaries: List[int] = []
    offset = 0

    for folder in folders:
        csv_path = folder / "rate.csv"
        if not csv_path.exists():
            continue

        df = pd.read_csv(csv_path)
        if df.empty or "second" not in df or "pair" not in df or "coincidences" not in df:
            continue

        # Ensure predictable ordering inside each run.
        df = df.sort_values(["second", "pair"]).reset_index(drop=True)
        df["abs_second"] = df["second"] + offset
        df["run"] = folder.name
        frames.append(df)

        last_sec = int(df["second"].max())
        offset += last_sec + 1  # +1 because seconds are inclusive
        boundaries.append(offset)

    combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return combined, boundaries


def compute_vis_qber(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute visibility and QBER per absolute second.
    Visibility/QBER (per second):
      HV basis: same = HH+VV, opp = HV+VH
      DA basis: same = DD+AA, opp = DA+AD
      vis = (same - opp) / (same + opp); qber = opp / (same + opp)
      Total = mean of HV and DA (ignoring missing basis).
    Returns wide DataFrame with vis/qber columns.
    """
    if df.empty:
        raise SystemExit("No rate.csv data found to plot.")

    pivot = df.pivot_table(index="abs_second", columns="pair", values="coincidences", aggfunc="sum").fillna(0)
    pivot = pivot.sort_index()

    def get(col):
        return pivot[col] if col in pivot else 0

    hh, vv, hv, vh = get("HH"), get("VV"), get("HV"), get("VH")
    dd, aa, da, ad = get("DD"), get("AA"), get("DA"), get("AD")

    total_coinc = pivot.sum(axis=1)

    # Drop low-count seconds from all downstream calculations
    mask_ok = total_coinc >= 40
    pivot = pivot[mask_ok]
    total_coinc = total_coinc[mask_ok]

    hh = pivot["HH"] if "HH" in pivot else 0
    vv = pivot["VV"] if "VV" in pivot else 0
    hv = pivot["HV"] if "HV" in pivot else 0
    vh = pivot["VH"] if "VH" in pivot else 0
    dd = pivot["DD"] if "DD" in pivot else 0
    aa = pivot["AA"] if "AA" in pivot else 0
    da = pivot["DA"] if "DA" in pivot else 0
    ad = pivot["AD"] if "AD" in pivot else 0

    hv_same = hh + vv
    hv_opp = hv + vh
    hv_total = hv_same + hv_opp
    vis_hv = (hv_same - hv_opp) / hv_total.replace(0, pd.NA)
    qber_hv = hv_opp / hv_total.replace(0, pd.NA)

    da_same = dd + aa
    da_opp = da + ad
    da_total = da_same + da_opp
    vis_da = (da_same - da_opp) / da_total.replace(0, pd.NA)
    qber_da = da_opp / da_total.replace(0, pd.NA)

    # Filter out negative visibility values (and their paired QBER)
    vis_hv = vis_hv.mask(vis_hv < 0)
    qber_hv = qber_hv.mask(vis_hv.isna())
    vis_da = vis_da.mask(vis_da < 0)
    qber_da = qber_da.mask(vis_da.isna())

    vis_total = pd.concat([vis_hv, vis_da], axis=1).mean(axis=1, skipna=True)
    qber_total = pd.concat([qber_hv, qber_da], axis=1).mean(axis=1, skipna=True)

    brightness = total_coinc / (0.8 * 1.58)  # given definition

    return pd.DataFrame(
        {
            "abs_second": pivot.index,
            "vis_hv": vis_hv,
            "qber_hv": qber_hv,
            "vis_da": vis_da,
            "qber_da": qber_da,
            "vis_total": vis_total,
            "qber_total": qber_total,
            "total_coinc": total_coinc,
            "brightness": brightness,
        }
    ).reset_index(drop=True)


def plot_vis_qber(metrics: pd.DataFrame, boundaries: List[int], out_path: Path) -> None:
    """Plot visibility and QBER (%) across absolute seconds."""
    if metrics.empty:
        raise SystemExit("No visibility/QBER data to plot.")

    m = metrics.copy()
    for col in m.columns:
        if col != "abs_second":
            m[col] = m[col].replace(pd.NA, np.nan).astype(float)

    def fmt_stats(series: pd.Series, percent=True) -> str:
        vals = series.dropna()
        if vals.empty:
            return "mean=nan, std=nan"
        mean = vals.mean()
        std = vals.std()
        if percent:
            mean *= 100
            std *= 100
        return f"mean={mean:.2f}%, std={std:.2f}%"

    t = m["abs_second"]
    fig, ax_vis = plt.subplots(figsize=(12, 6))
    ax_q = ax_vis.twinx()

    vis_lines = []
    q_lines = []

    vis_lines.append(
        ax_vis.plot(
            t,
            m["vis_total"] * 100,
            label=f"Total visibility ({fmt_stats(m['vis_total'])})",
            color="tab:blue",
        )[0]
    )
    q_lines.append(
        ax_q.plot(
            t,
            m["qber_total"] * 100,
            label=f"Total QBER ({fmt_stats(m['qber_total'])})",
            color="tab:red",
            linestyle="--",
        )[0]
    )

    # Optional basis-specific overlays
    if not m["vis_hv"].isna().all():
        vis_lines.append(
            ax_vis.plot(
                t,
                m["vis_hv"] * 100,
                label=f"HV visibility ({fmt_stats(m['vis_hv'])})",
                color="tab:green",
                alpha=0.5,
            )[0]
        )
    if not m["vis_da"].isna().all():
        vis_lines.append(
            ax_vis.plot(
                t,
                m["vis_da"] * 100,
                label=f"DA visibility ({fmt_stats(m['vis_da'])})",
                color="tab:orange",
                alpha=0.6,
            )[0]
        )
    if not m["qber_hv"].isna().all():
        q_lines.append(
            ax_q.plot(
                t,
                m["qber_hv"] * 100,
                label=f"HV QBER ({fmt_stats(m['qber_hv'])})",
                color="tab:brown",
                linestyle="--",
                alpha=0.5,
            )[0]
        )
    if not m["qber_da"].isna().all():
        q_lines.append(
            ax_q.plot(
                t,
                m["qber_da"] * 100,
                label=f"DA QBER ({fmt_stats(m['qber_da'])})",
                color="tab:pink",
                linestyle="--",
                alpha=0.6,
            )[0]
        )

    for boundary in boundaries[:-1]:
        ax_vis.axvline(boundary, color="gray", linestyle="--", alpha=0.25, linewidth=0.8)

    ax_vis.set_xlabel("Absolute second (concatenated runs)")
    ax_vis.set_ylabel("Visibility (%)")
    ax_q.set_ylabel("QBER (%)")
    ax_vis.set_title("Visibility and QBER across CoincEvents runs")
    ax_vis.grid(True, alpha=0.2)

    # Combined legend
    handles = vis_lines + q_lines
    labels = [h.get_label() for h in handles]
    ax_vis.legend(handles, labels, ncol=3, fontsize=8, loc="upper right")

    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    print(f"Wrote plot: {out_path}")


def plot_totals(metrics: pd.DataFrame, boundaries: List[int], out_path: Path) -> None:
    """Plot total coincidences and brightness on twin axes, with stats."""
    if metrics.empty:
        raise SystemExit("No totals data to plot.")

    m = metrics.copy()
    m["total_coinc"] = m["total_coinc"].astype(float)
    m["brightness"] = m["brightness"].astype(float)
    t = m["abs_second"]

    def fmt_stats(series: pd.Series) -> str:
        vals = series.dropna()
        if vals.empty:
            return "mean=nan, std=nan"
        return f"mean={vals.mean():.2f}, std={vals.std():.2f}"

    fig, ax_c = plt.subplots(figsize=(12, 6))
    ax_b = ax_c.twinx()

    line_c = ax_c.plot(t, m["total_coinc"], color="tab:blue",
                       label=f"Total coincidences ({fmt_stats(m['total_coinc'])})")[0]
    line_b = ax_b.plot(t, m["brightness"], color="tab:orange", linestyle="--",
                       label=f"Brightness ({fmt_stats(m['brightness'])})")[0]

    for boundary in boundaries[:-1]:
        ax_c.axvline(boundary, color="gray", linestyle="--", alpha=0.25, linewidth=0.8)

    ax_c.set_xlabel("Absolute second (concatenated runs)")
    ax_c.set_ylabel("Total coincidences")
    ax_b.set_ylabel("Brightness (cc/0.8/1.58)")
    ax_c.set_title("Total coincidences and brightness across CoincEvents runs")
    ax_c.grid(True, alpha=0.2)

    handles = [line_c, line_b]
    labels = [h.get_label() for h in handles]
    ax_c.legend(handles, labels, ncol=2, fontsize=8, loc="upper right")

    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    print(f"Wrote plot: {out_path}")


def main():
    # ----------------
    # Config
    # ----------------
    data_dir = Path("CoincEvents")           # directory containing run folders
    save_csv: Path | None = Path("combined_rate.csv")  # set to None to skip saving
    out_png = Path("vis_qber_consecutive.png")  # output plot
    out_totals_png = Path("totals_brightness_consecutive.png")  # totals/brightness plot
    plot_only_csv: Path | None = None        # if set, plot this CSV directly instead of aggregating
    plot_stride = 1                          # take every Nth point when plotting (>=1)

    if plot_only_csv is not None:
        csv_path = plot_only_csv
        if not csv_path.exists():
            raise SystemExit(f"CSV {csv_path} does not exist.")
        df = pd.read_csv(csv_path)
        if df.empty:
            raise SystemExit(f"CSV {csv_path} is empty.")
        if "abs_second" not in df.columns:
            if "second" not in df.columns:
                raise SystemExit("CSV must contain 'abs_second' or 'second' column.")
            df["abs_second"] = df["second"]
        boundaries: List[int] = []
    else:
        if not data_dir.exists():
            raise SystemExit(f"Data directory {data_dir} does not exist.")

        df, boundaries = load_rate_tables(data_dir)
        if df.empty:
            raise SystemExit(f"No rate.csv files found under {data_dir}")

        if save_csv is not None:
            df.to_csv(save_csv, index=False)
            print(f"Saved combined CSV: {save_csv}")

    metrics = compute_vis_qber(df)
    if plot_stride > 1:
        metrics = metrics.iloc[::plot_stride].reset_index(drop=True)
        boundaries = [b for b in boundaries if b % plot_stride == 0]
    plot_vis_qber(metrics, boundaries, out_png)
    plot_totals(metrics, boundaries, out_totals_png)


if __name__ == "__main__":
    main()
