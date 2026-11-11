# Release validation record (v1.0.0 candidate)

**Prepared:** 2026-11-06. This document records checks performed before any public GitHub/Zenodo/PyPI publication.

## Completed checks

- Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0, NetworkX 3.6.1, Matplotlib 3.10.8: **40 automated tests passed**, including the opt-in pixel-by-pixel comparison against the manuscript's original S4 figure on this specific renderer.
- Exact rational Markov certificate: both `8^4` four-block arrays have positive integer numerators and valid stationary three-block marginals; the ordinal-probability and lag-three-covariance equalities pass exact-integer checks.
- Recomputed Markov normalized OFS separation: `0.23122927994840325`; matching OPS and graph Fisher energy to floating-point precision.
- Analytical Gaussian AR(1): positive equal-energy parameter `0.6536965295755298`; Fisher-energy allocation TV distance `0.4127584219606162`.
- Executed both logistic-map and Markov simulation scripts with explicit reduced-size settings, confirming that they write output CSVs, metadata, figures, and consistent exact-reference quantities. **These reduced-size tests are not the full paper simulations.**
- Built a pure-Python wheel, verified the bundled rational certificate and fixed Cayley coordinates are included, and loaded the installed wheel from outside the source directory. The test environment reused preinstalled numerical dependencies rather than downloading them into an isolated venv.
- Built and inspected a source distribution containing the documentation, reproducibility scripts, tests, and public-release configuration.

- **External `ordpy` comparison (2026-11-11):** a Google Colab run reported `ordfisher` 1.0.0, `ordpy` 1.2.3, and Python 3.13.18; all 18 tested ordinal distributions matched exactly (maximum error 0), and `m=4` spectral consistency passed. This result is user-provided Colab output, not an independently captured machine-readable report. See `docs/ordpy-concordance.md`.

## Checks or decisions remaining before an actual public release

- The optional `ordpy` backend was not tested with the real dependency in the original local release environment; this gap was addressed by the external Colab concordance run recorded above. Additional independent test runs remain welcome.
- The full 510-realization logistic protocol and full-size 60-realization Markov simulation were run in earlier project packages, but **not independently rerun end-to-end in this v1.0.0 release candidate**. The available scripts preserve the documented full protocols.
- Only Python 3.13 was run in this container. The GitHub Actions workflow is configured to check Python 3.10, 3.11, 3.12 and 3.13 **after the repository becomes available**.
- `python -m twine check dist/*` is not recorded as completed in this environment (Twine is not installed); use the release workflow and checklist before publication.
- The four software authors and the maintainer have now been identified in the release metadata. All relevant rights holders still need to approve the final BSD-3-Clause license and copyright notice before publication.
- The GitHub repository address https://github.com/FranciscoX78/ordfisher is confirmed. A v1.0.0 tagged release, PyPI distribution, and Zenodo DOI remain to be verified; no DOI is invented.
- The exact Markov file is a **new reconstructed rational certificate**, not the unrecovered original certificate. All publication statements must preserve this distinction.
- Original historical simulation seeds were not recovered. Numerical reproductions use **new, recorded seeds**.

## Attribution and AI-disclosure update (2026-11-06)

- The four software authors are Giovanni Campanini, Fernando Pose, Carlos Garcia, and Francisco O. Redelico, as confirmed by the project lead; Francisco O. Redelico is the primary maintainer.
- Updated `pyproject.toml`, `CITATION.cff`, `AUTHORS.md`, the README and the BSD-3-Clause copyright line to reflect these names. Final approval of the license and attribution by all relevant rights holders remains necessary before public release.
- Added `docs/ai-use.md` to describe OpenAI ChatGPT (GPT-6 in the November 2026 revision stage), the software and writing tasks assisted, and the scientific verification responsibilities of the human authors.
- Re-ran the test suite with the strict image-comparison setting: **40 passed**.
- Built the wheel and source distribution in the available environment and verified the wheel's authors and maintainer metadata and a clean `--no-deps` installation/import from a separate directory. No public upload has been made.
- Publisher-specific methodology and declaration text is included in the separately supplied revised *Chaos* manuscript. Full historical versions of any earlier AI sessions should be confirmed by the authors prior to submission if required.
