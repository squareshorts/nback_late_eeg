# Secondary dynamical modeling

This directory contains the secondary modeling and model-discrimination artifacts used to characterize the already-established late EEG power reduction. The primary whole-domain inference remains in the repository root and `scripts/`.

## Retained manuscript model

The retained model is a three-parameter delayed critically damped transient:

```text
g(t) = 0, t < delta
g(t) = -A * u * exp(1-u), u=(t-delta)/tau, t >= delta
```

It is used as a **phenomenological amplitude/gain description**, not as evidence for a specific cellular circuit, Hopf bifurcation, or oscillator topology.

Deterministic refit to `results/group_late_contrast_trajectory.csv` gives:

- A = 1.159 dB
- delay = 481.3 ms
- tau = 114.6 ms
- peak time = 595.9 ms
- group R2 = 0.881
- group RMSE = 0.126 dB
- LOSO median predictive R2 vs zero = 0.623
- positive predictive R2 in 13/16 participants
- LOSO median correlation = 0.659
- LOSO median RMSE = 0.500 dB

The unrestricted training-mean trajectory has median held-out RMSE 0.497 dB; the one-parameter constant model has median RMSE 0.538 dB.

## Phase-control result

Direct phase controls did not show a comparably robust higher-demand reduction:

- code-balanced PPC over 12 frozen max-statistic voxels: -0.0030, 95% CI [-0.0267, 0.0206], exact p=.798
- PPC at the F3/15-Hz core: -0.0251, 95% CI [-0.0571, 0.0069], exact p=.113
- spatial phase order R at 15 Hz, 620--720 ms: -0.00374, 95% CI [-0.00961, 0.00214], exact p=.191

These controls motivated describing the reproducible macroscopic effect primarily as amplitude/gain dynamics rather than phase-coherence loss.

## Reproduction

From the repository root:

```bash
python modeling/scripts/07_delayed_gain_model.py
```

This reads only committed derived CSVs and regenerates the titleless model figure.

No raw EEG or identifying participant data are included here.
