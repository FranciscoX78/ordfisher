# Changelog

## 1.0.0 (2026-11-06)

First public-release candidate of `ordfisher`.

- Identified the four named software authors and Francisco O. Redelico as maintainer.
- Documented substantive OpenAI ChatGPT assistance and human verification responsibilities.

- Calculate ordinal pattern distributions and the Ordinal Power Spectrum (OPS) and Ordinal Fisher Spectrum (OFS) on the adjacent-rank Cayley graph.
- Compute normalized OPS/OFS profiles, Cayley Fisher energy, spectral Fisher-energy allocation, and summary statistics.
- Provide an independent NumPy counting backend and an optional `ordpy` adapter for tie-free series.
- Build and render Cayley graphs, including the fixed layout matching the manuscript's S4 appendix figure.
- Include exact rational verification data for the reconstructed stationary Markov examples.
- Include deterministic AR(1) calculations and seeded logistic-map/Markov reproduction scripts.
- Add release metadata, BSD-3-Clause license, English documentation, continuous-integration configuration, and tests.

**Provenance:** the bundled rational certificate is a newly reconstructed exact realization of the paper's stated constraints, **not** the original historical certificate. Logistic and Markov simulations use recorded new seeds, not the original unrecovered simulation streams.
