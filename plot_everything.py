import numpy as np, pandas as pd, sys
import matplotlib.pyplot as plt
from pathlib import Path

# ===================================
# 1. Configuration
# ===================================
Setup_4 = False  # Change to True for 4-detector setup

data_dir = Path("Delay_Scan_Data")
names = ["delay_ns", "coinc"]
end_time = sys.argv[1]
seconds = range(0, int(end_time))

# ===================================
# 2. Detector Pair Setup
# ===================================
if Setup_4:
    # 4-detector setup: H1, V1, V2, H2
    pairs = {
        "H1–H2": (1, 4),
        "V1–V2": (2, 3),
    }
    same_HV = [(1, 4), (2, 3)]
    opp_HV  = [(1, 3), (2, 4)]
    same_DA, opp_DA = [], []  # no D/A detectors in this setup
else:
    # 8-detector setup
    pairs = {
        "H1–H2": (1, 5),
        "V1–V2": (2, 6),
        "D1–D2": (3, 7),
        "A1–A2": (4, 8),
    }
    same_HV = [(1,5),(2,6)]
    opp_HV  = [(1,6),(2,5)]
    same_DA = [(3,7),(4,8)]
    opp_DA  = [(3,8),(4,7)]

# ===================================
# 3. Helper Functions
# ===================================
def get_peak_delay_and_count(i, j, k):
    """Return (delay_at_max, max_count) for a given pair and second."""
    f = data_dir / f"delay_scan_{i}_vs_{j}_second_{k}.csv"
    if not f.exists():
        return None, 0
    df = pd.read_csv(f, names=names)
    if df.empty:
        return None, 0
    idx_max = df["coinc"].idxmax()
    delay_at_max = df.loc[idx_max, "delay_ns"]
    max_count = df.loc[idx_max, "coinc"]
    return delay_at_max, max_count


def get_count_at_delay(i, j, k, delay_target):
    """Return coincidence count at a given delay (nearest value)."""
    f = data_dir / f"delay_scan_{i}_vs_{j}_second_{k}.csv"
    if not f.exists():
        return 0
    df = pd.read_csv(f, names=names)
    if df.empty:
        return 0
    # Find the coincidence closest to the target delay
    idx = (df["delay_ns"] - delay_target).abs().idxmin()
    return df.loc[idx, "coinc"]


def counts_same_opp(k, same_pairs, opp_pairs):
    """
    For each same/opp pair, ensure opposite coincidence is evaluated
    at the delay where the same pair peaks.
    """
    C_same_total, C_opp_total = 0, 0
    for (same_i, same_j), (opp_i, opp_j) in zip(same_pairs, opp_pairs):
        delay_peak, C_same = get_peak_delay_and_count(same_i, same_j, k)
        if delay_peak is None:
            continue
        C_opp = get_count_at_delay(opp_i, opp_j, k, delay_peak)
        C_same_total += C_same
        C_opp_total += C_opp
    return C_same_total, C_opp_total


def compute_metrics(same_pairs, opp_pairs):
    """Compute visibility, contrast, QBER per second."""
    V, CR, Q, C_same_all, C_opp_all = [], [], [], [], []
    for k in seconds:
        C_same, C_opp = counts_same_opp(k, same_pairs, opp_pairs)
        C_same_all.append(C_same)
        C_opp_all.append(C_opp)
        T = C_same + C_opp
        if T == 0 or C_opp == 0:
            V.append(np.nan); CR.append(np.nan); Q.append(np.nan)
            continue
        v  = (C_same - C_opp) / T
        cr = C_same / C_opp
        V.append(v)
        CR.append(cr)
        Q.append(C_opp / T)
    return np.array(V), np.array(CR), np.array(Q), np.array(C_same_all), np.array(C_opp_all)

# ===================================
# 4. Plot All Delay Scans
# ===================================
global_max = 0
for (i, j) in pairs.values():
    for k in seconds:
        f = data_dir / f"delay_scan_{i}_vs_{j}_second_{k}.csv"
        if not f.exists():
            continue
        df = pd.read_csv(f, names=["delay_ns", "coincidences"])
        if not df.empty:
            global_max = max(global_max, df["coincidences"].max())

fig, axes = plt.subplots(len(pairs), 1, figsize=(8, 6), sharex=True)
if len(pairs) == 1:
    axes = [axes]
for ax, (label, (i, j)) in zip(axes, pairs.items()):
    for k in seconds:
        f = data_dir / f"delay_scan_{i}_vs_{j}_second_{k}.csv"
        if not f.exists():
            continue
        df = pd.read_csv(f, names=["delay_ns", "coincidences"])
        ax.plot(df["delay_ns"], df["coincidences"], label=f"s {k+1}")
    ax.set_ylim(0, global_max * 1.1)
    ax.set_title(label)
    ax.set_ylabel("Coincidences")
    ax.legend(fontsize=8)
    ax.grid(True)
axes[-1].set_xlabel("Relative delay (ns)")
plt.tight_layout()
plt.show()

# ===================================
# 5. Compute Metrics
# ===================================
V_HV, CR_HV, Q_HV, Cs_HV, Co_HV = compute_metrics(same_HV, opp_HV)

if Setup_4:
    V_DA = CR_DA = Q_DA = np.full_like(V_HV, np.nan)
else:
    V_DA, CR_DA, Q_DA, Cs_DA, Co_DA = compute_metrics(same_DA, opp_DA)

t = list(seconds)
individual_coinc = {}

# ===================================
# 6. Individual Coincidence Trends
# ===================================
if Setup_4:
    fig, ax = plt.subplots(figsize=(8, 5))
    combos = {
        "H1–H2 (same)": (1, 4),
        "V1–V2 (same)": (2, 3),
        "H1–V2 (opp)":  (1, 3),
        "V1–H2 (opp)":  (2, 4),
    }
    for label, (i, j) in combos.items():
        coinc_series = [get_peak_delay_and_count(i, j, k)[1] for k in seconds]
        individual_coinc[label] = coinc_series
        style = 'o-' if "same" in label else 'o--'
        ax.plot(t, coinc_series, style, label=label)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Max coincidences")
    ax.set_title("H/V Basis Coincidences vs Time")
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.show()

else:
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharex=True, sharey=True)

    combos_HV = {
        "H1–H2 (same)": (1, 5),
        "V1–V2 (same)": (2, 6),
        "H1–V2 (opp)":  (1, 6),
        "V1–H2 (opp)":  (2, 5),
    }
    combos_DA = {
        "D1–D2 (same)": (3, 7),
        "A1–A2 (same)": (4, 8),
        "D1–A2 (opp)":  (3, 8),
        "A1–D2 (opp)":  (4, 7),
    }

    for ax, (title, combos) in zip(axes, [("H/V Basis", combos_HV), ("D/A Basis", combos_DA)]):
        for label, (i, j) in combos.items():
            coinc_series = [get_peak_delay_and_count(i, j, k)[1] for k in seconds]
            individual_coinc[label] = coinc_series
            style = 'o-' if "same" in label else 'o--'
            ax.plot(t, coinc_series, style, label=label)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Max coincidences")
        ax.set_title(f"{title} Coincidences vs Time")
        ax.legend()
        ax.grid(True)

    plt.tight_layout()
    plt.show()

# ===================================
# 7. Visibility, Contrast, QBER Plots
# ===================================
fig, ax = plt.subplots(2, 2, figsize=(9, 8), sharex=True)
axes = ax.ravel()

axes[0].plot(t, V_HV * 100, 'o-', label="HV", color='tab:blue')
if not Setup_4:
    axes[0].plot(t, V_DA * 100, 'o--', label="DA", color='tab:orange')
axes[0].set_ylabel("Visibility (%)")
axes[0].legend()

axes[1].plot(t, CR_HV, 'o-', label="HV", color='tab:blue')
if not Setup_4:
    axes[1].plot(t, CR_DA, 'o--', label="DA", color='tab:orange')
axes[1].set_ylabel("Contrast ratio")

axes[2].plot(t, Q_HV * 100, 'o-', label="HV", color='tab:blue')
if not Setup_4:
    axes[2].plot(t, Q_DA * 100, 'o--', label="DA", color='tab:orange')
axes[2].set_ylabel("QBER (%)")
axes[2].set_xlabel("Time (s)")
axes[2].legend()

axes[3].axis("off")
for a in axes:
    a.grid(True)

plt.suptitle("Visibility, Contrast, and QBER vs Time", y=1.02)
plt.tight_layout()
plt.show()

# ===================================
# 8. Summary Statistics
# ===================================
summary_data = {
    "HV_Vis_mean": np.nanmean(V_HV), "HV_Vis_std": np.nanstd(V_HV),
    "HV_QBER_mean": np.nanmean(Q_HV), "HV_QBER_std": np.nanstd(Q_HV),
    "HV_Contrast_mean": np.nanmean(CR_HV), "HV_Contrast_std": np.nanstd(CR_HV),
}

if not Setup_4:
    summary_data.update({
        "DA_Vis_mean": np.nanmean(V_DA), "DA_Vis_std": np.nanstd(V_DA),
        "DA_QBER_mean": np.nanmean(Q_DA), "DA_QBER_std": np.nanstd(Q_DA),
        "DA_Contrast_mean": np.nanmean(CR_DA), "DA_Contrast_std": np.nanstd(CR_DA),
    })

# Add per-pair coincidence statistics
for label, coincs in individual_coinc.items():
    summary_data[f"{label}_mean"] = np.nanmean(coincs)
    summary_data[f"{label}_std"] = np.nanstd(coincs)

summary = pd.Series(summary_data).to_frame("Value")
# print(f"HV QBER: {Q_HV}")
# print(f"DA QBER: {Q_DA}")

print("\n--- Summary Statistics ---")
print(summary.round(4))
