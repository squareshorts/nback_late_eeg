#!/usr/bin/env python3
"""Regenerate the titleless supplementary phase-control figure from safe derived results."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
FIGURES = HERE.parent / "figures"


def main():
    df = pd.read_csv(RESULTS / "phase_power_control_summary.csv")
    phase = df[df["analysis"].isin(["Phase PPC (code-balanced)", "Spatial phase order R"])].copy()
    labels = [
        "PPC\n12 max-stat voxels",
        "PPC\nF3/15-Hz core",
        "Spatial phase order\n15 Hz, 620--720 ms",
    ]
    x = np.arange(len(phase))
    y = phase["estimate"].to_numpy(float)
    yerr = np.vstack([
        y - phase["ci95_low"].to_numpy(float),
        phase["ci95_high"].to_numpy(float) - y,
    ])

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.errorbar(x, y, yerr=yerr, fmt="o", capsize=4)
    ax.axhline(0, linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Higher-demand minus 0-back phase metric")
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    FIGURES.mkdir(exist_ok=True)
    fig.savefig(FIGURES / "FigureS3_phase_controls.svg", bbox_inches="tight")


if __name__ == "__main__":
    main()
