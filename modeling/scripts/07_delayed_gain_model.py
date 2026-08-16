#!/usr/bin/env python3
"""Reproduce the delayed critically damped amplitude-gating analysis.

Inputs are safe derived CSVs in ../results. No raw EEG is read.
The model is phenomenological and is not a cellular or bifurcation model.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
FIGURES = HERE.parent / "figures"

def response(t, A, delta, tau):
    t = np.asarray(t, dtype=float)
    y = np.zeros_like(t)
    mask = t >= delta
    u = (t[mask] - delta) / tau
    y[mask] = -A * u * np.exp(1.0 - u)
    return y

def main():
    df = pd.read_csv(RESULTS / "group_late_contrast_trajectory.csv")
    t = df["time_ms"].to_numpy(float)
    y = df["empirical_mean_step_db"].to_numpy(float)
    popt, _ = curve_fit(response, t, y, p0=[1.15, 480.0, 115.0], bounds=([0.0, 350.0, 20.0], [5.0, 600.0, 500.0]), maxfev=100000)
    fit = response(t, *popt)
    rmse = np.sqrt(np.mean((y - fit) ** 2))
    r2 = 1.0 - np.sum((y - fit) ** 2) / np.sum((y - y.mean()) ** 2)
    print(f"A={popt[0]:.6f} dB")
    print(f"delta={popt[1]:.6f} ms")
    print(f"tau={popt[2]:.6f} ms")
    print(f"peak={popt[1]+popt[2]:.6f} ms")
    print(f"R2={r2:.6f}")
    print(f"RMSE={rmse:.6f} dB")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(t, y, marker="o", label="Empirical mean contrast")
    if "empirical_sem_db" in df:
        sem = df["empirical_sem_db"].to_numpy(float)
        ax.fill_between(t, y-sem, y+sem, alpha=0.2)
    ax.plot(t, fit, linewidth=2, label="Delayed critically damped model")
    ax.axhline(0, linewidth=1)
    ax.set_xlabel("Time from stimulus onset (ms)")
    ax.set_ylabel("Higher-demand minus 0-back power (dB)")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    FIGURES.mkdir(exist_ok=True)
    fig.savefig(FIGURES / "Figure4_dynamical_gain_model.svg", bbox_inches="tight")

if __name__ == "__main__":
    main()
