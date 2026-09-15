# N-back late EEG analysis: reproducibility archive

This repository is the frozen reproducibility archive for the final robustness analysis underlying the manuscript **“Distributed late stimulus-related spectral-power reduction accompanies higher-demand visuospatial N-back performance.”**

The archive contains the deterministic stimulus-balanced whole-domain analysis, leave-one-subject-out (LOSO) extraction, nontarget/response-censoring sensitivity analyses, integrity checks, and the retained secondary dynamical/phase-control analyses. Raw identifiable EEG recordings are not included. The committed safe intermediate data are sufficient to rerun the archived inferential pipeline without the legacy development directories.

## Final validation anchors

A clean reproduction of the frozen pipeline recovered the final manuscript-level results:

- 100,000-permutation whole-domain inference completed successfully.
- The deterministic stimulus-balanced late component survived family-wise correction; exactly **12 voxels** also survived voxelwise maximum-statistic correction.
- Correct-rejection nontarget LOSO late contrast: **-0.811 dB**, 95% CI **[-1.146, -0.477]**, **dz = -1.29**.
- All technically valid nontargets: **-0.818 dB**, 95% CI **[-1.149, -0.487]**.
- Response-censored all-nontarget analysis: **-0.812 dB**, 95% CI **[-1.143, -0.481]**.

The primary deterministic-balanced cluster result is stored in `DETERMINISTIC_STIMULUS_BALANCED_CLUSTER_RESULTS.csv`; the max-statistic result is stored in `DETERMINISTIC_STIMULUS_BALANCED_MAXSTAT.csv`; LOSO and trial-population sensitivity results are stored in the corresponding CSV files at repository root and under `data/`.

## Repository scope

This repository is an **analysis/reproducibility archive**, not the manuscript source repository. The `figures/` directory contains figures produced by the archived robustness/reproduction pipeline. Publication-layout figures that were assembled separately for the manuscript are not implied to have their original plotting scripts here unless an explicit regeneration script is committed.

The purpose of this boundary is to keep the public archive minimal, deterministic, and free of legacy working files.

## Environment

Python dependencies are listed in:

```text
environment/requirements.txt
```

Using `uv`:

```bash
uv venv .venv
uv pip install --python .venv -r environment/requirements.txt
```

Using standard Python/pip:

```bash
python -m venv .venv
python -m pip install -r environment/requirements.txt
```

Activate the environment in the usual way for your platform, or call the environment's Python executable directly.

## Primary execution order

Run from the repository root:

```bash
python scripts/01_build_balanced_maps.py
python scripts/02_balanced_whole_domain.py
python scripts/03_balanced_loso.py
python scripts/04_response_and_selection.py
python scripts/05_synthesize_final.py
python scripts/06_integrity.py
```

`02_balanced_whole_domain.py` performs the permutation-based whole-domain inference and is the computationally intensive stage.

The pipeline uses committed safe intermediate data rather than raw identifiable EEG files, then reruns the permutation inference, LOSO statistics, sensitivity analyses, synthesis, and integrity checks represented in this frozen archive.

## Secondary dynamical modeling and phase controls

The manuscript's secondary temporal model and phase-control artifacts are stored in `modeling/`. These analyses do not alter the primary whole-domain inference.

To reproduce the retained delayed critically damped temporal fit from the committed derived trajectory:

```bash
python modeling/scripts/07_delayed_gain_model.py
```

The retained group fit was summarized by **A = 1.159 dB**, **delay = 481.3 ms**, **time constant = 114.6 ms**, and **R2 = .881**. Leave-one-subject-out prediction yielded a median predictive **R2 = .623**, with positive predictive R2 in **13/16** participants.

See `modeling/README.md` and `modeling/MODEL_SELECTION.md` for the retained model, LOSO validation, phase controls, and documentation of exploratory dynamical analyses that were not retained as manuscript mechanisms.

## Integrity and interpretation

- `HASH_MANIFEST.csv` records the archive integrity manifest.
- `FINAL_ROBUSTNESS_EXECUTIVE_SUMMARY.md` records the final robustness interpretation and important scope limitations.
- The task blocks were administered in a fixed 0-back -> 1-back -> 2-back sequence; the archived result is therefore interpreted as a condition/sequence-associated late spectral-power difference rather than a causal estimate of parametric working-memory load.
- The early effect is not treated as a successful family-wise-confirmed replication under the final deterministic-balanced primary analysis.

## Reuse

When reusing this archive, preserve the deterministic stimulus-balancing procedure, the frozen whole-domain inference settings, and the distinction between primary inference and secondary descriptive modeling. Do not substitute legacy development outputs for the committed final files without documenting the change.
