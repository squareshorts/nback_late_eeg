# Model selection record

The modeling work deliberately separated **generative plausibility** from **empirical support**. Only the delayed critically damped amplitude/gain model is retained in the manuscript.

## 1. Hopf/Stuart-Landau crossing: rejected

A corrected stochastic Stuart-Landau network was tested with proper Euler-Maruyama noise scaling, an actual small-world topology, baseline equilibration, repeated realizations, and a post-stimulus sweep of the bifurcation parameter.

Key result: with baseline `a=+0.5`, the originally proposed crossing to `a=-0.5` produced a late contrast around -3.37 dB, far larger than the empirical -0.811 dB. The empirical magnitude was instead reproduced near `a_post=+0.20` to `+0.25`, while the oscillator remained on the positive-a limit-cycle side. Therefore the observed power effect does not identify or require a Hopf crossing.

## 2. Isochron/shear mechanism: not retained

Curved-isochron simulations showed that phase-amplitude shear can in principle convert transient radial perturbations into phase dispersion, but broad predeclared parameter sweeps were not predictive enough: the median late effect was about -0.17 dB, and only a minority of parameter sets fell near the empirical magnitude.

More importantly, direct empirical phase controls at locations fixed independently by the power result did not show a reliable higher-demand phase-coherence decrease (PPC or spatial phase order). This makes an isochron-driven phase-dispersion explanation speculative for the present scalp EEG.

## 3. Delayed amplitude/gain transient: retained

The retained model was selected because it describes the actual temporal morphology, uses only three interpretable parameters, and generalizes to held-out participants. It is intentionally phenomenological and should not be interpreted as a unique cellular mechanism.

Exploratory numerical files are preserved under `exploratory/` for provenance only and are **not manuscript evidence**.
