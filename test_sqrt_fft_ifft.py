"""
Smoke test for `sqrt_fft_ifft` using a synthetic waveform.

We build a histogram-like dataset whose bin populations follow a sum of
sinusoids/cosines at 20, 40, 80, and 100 Hz (cycles over the bin index).
Feeding those samples into `sqrt_fft_ifft` should yield a reconstructed
waveform whose spectrum contains the same low-frequency peaks (phase is
discarded inside the function, so we only compare magnitudes).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# Use non-interactive backend to avoid GUI requirements in tests.
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Prefer the build/ copy (contains sqrt_fft_ifft) if it exists; fall back to repo root.
repo_root = Path(__file__).resolve().parent
build_dir = repo_root / "build"
if build_dir.exists():
    sys.path.insert(0, str(build_dir))

from Nicely_Plotted_various_diffs import sqrt_fft_ifft  # noqa: E402


def generate_values(n_bins: int = 2048) -> np.ndarray:
    """Return repeated bin-center samples that encode target frequencies."""
    # Frequencies (cycles across n_bins) and relative amplitudes chosen to make
    # the peaks visible after the magnitude-only pipeline in sqrt_fft_ifft.
    components = [(20, 2.0), (40, 1.0), (80, 1.0), (100, 2.0)]
    idx = np.arange(n_bins, dtype=float)
    waveform = sum(
        amp * (np.sin(2 * np.pi * freq * idx / n_bins) + np.cos(2 * np.pi * freq * idx / n_bins))
        for freq, amp in components
    )
    # Shift to all-positive so the counts are valid histogram populations.
    counts = np.round((waveform - waveform.min() + 1.0) * 50).astype(int)
    # Bin centers on a uniform grid; sqrt_fft_ifft will recompute similar centers internally.
    centers = np.linspace(0.0, 1.0, n_bins, endpoint=False)
    return np.repeat(centers, counts)


def find_peaks(freqs: np.ndarray, magnitudes: np.ndarray, k: int = 6, fmax: float = 200.0):
    """Return the k largest peak frequencies below fmax, sorted ascending."""
    mask = freqs <= fmax
    trimmed_mag = magnitudes.copy()
    trimmed_mag[~mask] = 0
    trimmed_mag[0] = 0  # drop DC
    peak_idx = trimmed_mag.argsort()[-k:][::-1]
    return np.sort(freqs[peak_idx])


def test_sqrt_fft_ifft_recovers_frequencies():
    values = generate_values()
    centers, f_rec = sqrt_fft_ifft(values, bins=2048)
    assert centers is not None and f_rec is not None, "Function returned None outputs"

    dt = centers[1] - centers[0]
    freqs = np.fft.rfftfreq(len(f_rec), d=dt)
    mag = np.abs(np.fft.rfft(f_rec))

    peaks = find_peaks(freqs, mag, k=8, fmax=200.0)
    # Expect the target frequencies to appear within a small tolerance.
    expected = np.array([20.0, 40.0, 80.0, 100.0])
    for f in expected:
        assert np.any(np.isclose(peaks, f, atol=1.5)), f"Missing peak near {f} Hz (found {peaks})"

    plot_diagnostics(values, centers, f_rec, freqs, mag)


def plot_diagnostics(values: np.ndarray, centers, f_rec, freqs, mag):
    """Save quick-look plots: input histogram, reconstructed signal, and FFT magnitude."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    # Original histogram proxy (values are repeated bin centers).
    axes[0].hist(values, bins=200, color="steelblue", alpha=0.7)
    axes[0].set_title("Original Samples (values)")
    axes[0].set_xlabel("Sample")
    axes[0].set_ylabel("Count")

    # Reconstructed time-domain waveform.
    axes[1].plot(centers, f_rec, color="darkorange", lw=1.2)
    axes[1].set_title("sqrt_fft_ifft Output")
    axes[1].set_xlabel("Time bin")
    axes[1].set_ylabel("Amplitude (arb.)")

    # FFT magnitude of reconstructed waveform.
    axes[2].plot(freqs, mag, color="seagreen", lw=1.0)
    axes[2].set_xlim(0, 150)
    axes[2].set_title("FFT Magnitude of Reconstruction")
    axes[2].set_xlabel("Frequency (Hz)")
    axes[2].set_ylabel("Magnitude")

    fig.tight_layout()
    # Write inside current working directory to stay sandbox-friendly.
    out_path = Path.cwd() / "sqrt_fft_ifft_diagnostics.png"
    fig.savefig(out_path, dpi=150)


if __name__ == "__main__":
    test_sqrt_fft_ifft_recovers_frequencies()
    print("sqrt_fft_ifft frequency smoke test: PASS")
