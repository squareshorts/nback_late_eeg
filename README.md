# N-back Late EEG Analysis Reproducibility

This repository contains the exact frozen analysis and scripts for the n-back late cluster results.

## Execution Order
1. 01_build_balanced_maps.py
2. 02_balanced_whole_domain.py
3. 03_balanced_loso.py
4. 04_response_and_selection.py
5. 05_synthesize_final.py
6. 06_integrity.py

## Running the Analysis
Activate the virtual environment inside `environment/` and run the scripts in order. The pipeline bypasses raw identifiable EEG extraction by loading the intermediate safe data, then reruns the permutation-based inferences and final statistics exactly as reported in the final robustness pass.

## Secondary dynamical modeling

The manuscript's secondary temporal model and phase-control artifacts are stored in `modeling/`. This layer does not change the primary whole-domain inference. It contains only safe derived results, model-selection documentation, figure-regeneration code, and a deterministic refit of the three-parameter delayed critically damped amplitude/gain model.

To reproduce the retained dynamical fit from the committed derived trajectory:

```bash
python modeling/scripts/07_delayed_gain_model.py
```

See `modeling/README.md` and `modeling/MODEL_SELECTION.md` for the retained model, LOSO validation, phase controls, and the exploratory Hopf/isochron analyses that were not retained as manuscript mechanisms.
