# v1.0.0 — Frozen reproducibility archive

This release freezes the reproducibility archive for the final robustness analysis underlying the manuscript **“Distributed late stimulus-related spectral-power reduction accompanies higher-demand visuospatial N-back performance.”**

## Reproduced final results

The archived pipeline was reproduced in a clean, self-contained environment without requiring the legacy `nback_TEMP` development tree.

- 100,000-permutation whole-domain inference completed successfully.
- Exactly **12 voxels** survived voxelwise maximum-statistic family-wise correction within the retained late component.
- Correct-rejection nontarget LOSO late contrast: **-0.811 dB**, 95% CI **[-1.146, -0.477]**, **dz = -1.29**.
- All technically valid nontargets: **-0.818 dB**, 95% CI **[-1.149, -0.487]**.
- Response-censored all-nontarget analysis: **-0.812 dB**, 95% CI **[-1.143, -0.481]**.

## Archive contents

The release contains:

- deterministic equal-weight stimulus-balanced whole-domain inference;
- 100,000-permutation cluster and max-statistic correction;
- leave-one-subject-out feature extraction and participant-level contrasts;
- all-nontarget and response-censored sensitivity analyses;
- integrity/hash checks;
- retained secondary delayed critically damped temporal modeling;
- phase-control analyses and their reproducibility artifacts;
- safe derived/intermediate data sufficient to rerun the frozen analysis without raw identifiable EEG recordings.

## Scope

This repository is the frozen **analysis/reproducibility archive**. It is not the manuscript-source repository, and publication-layout figures assembled elsewhere are not assumed to have their original plotting scripts here unless explicitly committed.

The fixed 0-back -> 1-back -> 2-back block order remains an experimental limitation. The retained late effect should therefore be interpreted as condition/sequence-associated rather than as a causal estimate of parametric working-memory load.

The earlier broad early-effect claim is not retained as a family-wise-confirmed result under the final deterministic-balanced primary analysis.

## Reproduction

From the repository root, install the dependencies in `environment/requirements.txt` and run:

```bash
python scripts/01_build_balanced_maps.py
python scripts/02_balanced_whole_domain.py
python scripts/03_balanced_loso.py
python scripts/04_response_and_selection.py
python scripts/05_synthesize_final.py
python scripts/06_integrity.py
```

For the retained secondary temporal model:

```bash
python modeling/scripts/07_delayed_gain_model.py
```

See `README.md`, `FINAL_ROBUSTNESS_EXECUTIVE_SUMMARY.md`, `modeling/README.md`, and `modeling/MODEL_SELECTION.md` for details.
