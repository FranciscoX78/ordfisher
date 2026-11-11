# ordfisher

**Ordinal Power Spectrum (OPS) and Ordinal Fisher Spectrum (OFS) in Python**

`ordfisher` provides reproducible, research-oriented implementations of two graph spectra of ordinal-pattern probability distributions. Both use the adjacent-rank Cayley graph of the symmetric group; OPS analyzes the probability vector, while OFS analyzes its componentwise square root. The library also supplies the associated graph-based Fisher energy, spectrum summaries, and publication-quality Cayley-graph drawings.

> **Release status:** source prepared for the public repository [FranciscoX78/ordfisher](https://github.com/FranciscoX78/ordfisher). A GitHub **v1.0.0 release**, PyPI publication, and a Zenodo DOI must not be claimed until each has been independently verified. See [the release checklist](docs/release-checklist.md).

## Features

- Compute all `m!` ordinal-pattern probabilities (including missing patterns), with explicit delay and tie policies.
- Build the **adjacent-rank** Cayley graph and calculate its complete normalized-Laplacian spectrum.
- Obtain OPS and OFS powers grouped by the entire eigenspace of each distinct eigenvalue.
- Compute scalar Cayley Fisher energy, normalized OFS/OPS profiles, Fisher-energy allocation, and spectral centroids.
- Visualize the Cayley graph. For `m=4`, the default `appendix` layout reproduces the layout and styling of the OFS manuscript's S4 figure.
- Verify spectral identities and reconstruct analytical Gaussian AR(1) results.
- Reproduce the logistic-map and stationary Markov examples with recorded random seeds; validate a rational Markov certificate using exact integer arithmetic.
- Optionally use `ordpy` to produce ordinal-pattern probabilities **when the input contains no ties**.

## Installation

Requires **Python 3.10 or newer**. The package uses NumPy, SciPy, Matplotlib, and NetworkX. Install from the downloaded repository or source archive:

```bash
python -m pip install .
```

For development and tests:

```bash
python -m pip install -e '.[test]'
python -m pytest -q
```

Optional dependencies can be installed separately:

```bash
python -m pip install '.[ordpy]'    # independent ordinal-counting adapter
python -m pip install '.[fast]'     # optional Numba acceleration
python -m pip install '.[rebuild]'  # optional rational-certificate reconstruction
```

The library has **not been verified as available through PyPI**, so `pip install ordfisher` is not a verified release instruction at this stage. See the [release checklist](docs/release-checklist.md) before publishing.

## Quick start: analyze a time series

```python
import numpy as np
from ordfisher import spectra_from_series, check_spectral_identities

x = np.random.default_rng(123).normal(size=10_000)
result = spectra_from_series(x, m=4, tau=1)

print("Graph frequencies:", result.frequencies)
print("OPS powers:", result.ops)
print("OFS powers:", result.ofs)
print("Cayley Fisher energy:", result.fisher_energy)
print("OPS-OFS profile distance:", result.ops_ofs_tv())

# Verify independent mathematical identities against this probability vector.
print(check_spectral_identities(result.probabilities, m=4))
```

**Interpret the three normalized profiles carefully:**

- `result.ops_profile` distributes OPS power over **positive** graph frequencies.
- `result.ofs_profile` distributes OFS power over **positive** graph frequencies.
- `result.fisher_profile` distributes **Fisher energy** over graph frequencies, weighting OFS by the eigenvalue.

The arrays `result.ops` and `result.ofs` *include* the zero-frequency term. Normalized profiles are `None` when their corresponding positive-frequency power is numerically zero. Frequencies are **graph frequencies, not temporal Fourier frequencies**.

## Start from a known ordinal distribution

```python
from ordfisher import spectra_from_probabilities

# All 3! probabilities, in lexicographic order of argsort permutations:
p = [1/3, 1/15, 2/15, 1/5, 4/15, 0]
r = spectra_from_probabilities(p, m=3)
print(r.ofs, r.ops, r.fisher_energy)
```

If you supply the probabilities directly, `p` **must** contain all `m!` patterns, including zero entries, in **lexicographic order** of time-index permutations sorted by their observed values (`argsort`). The function checks normalization and does not silently renormalize invalid probabilities. For input series, `ties='stable'` breaks exact ties using the earlier time index; `ties='raise'` rejects tied windows. The optional `engine='ordpy'` adapter is restricted to **tie-free** data to avoid conflicting tie conventions.

## Draw the Cayley graph

```python
from ordfisher import cayley_graph, plot_cayley_graph

G = cayley_graph(4)
assert (G.number_of_nodes(), G.number_of_edges()) == (24, 36)

fig, ax = plot_cayley_graph(4, layout="appendix", save="cayley_S4.pdf")
# For another layout: layout="spring" or "kamada_kawai".
```

For `m=4`, `plot_cayley_graph(4)` uses the fixed appendix layout by default. The layout, edges, and colors are fixed; pixel-perfect equivalence across operating systems is not guaranteed because fonts and Matplotlib renderers differ. The **same mathematical graph** is used in the spectral calculations.

## Supported scope

| Operation | Supported pattern length |
|---|---|
| Full OPS/OFS (dense complete spectrum) | `2 <= m <= 6` |
| Build/draw Cayley graph without spectral decomposition | `2 <= m <= 7` |
| Analytical stationary Gaussian AR(1) ordinal distribution | `m = 4` |
| Reconstructed rational Markov example | `m = 4` |

Full spectral decomposition is dense and requires `m!`-by-`m!` matrices. The package deliberately **does not** return truncated eigenvalues as if they represented the complete spectrum. For large `m`, time and memory requirements increase rapidly.

## Reproduce the companion paper

The scripts are separated from the reusable numerical library. From the repository root:

```bash
python scripts/verify_markov.py --out results/markov/verification.json
python scripts/verify_ar1.py --out results/ar1
python scripts/run_logistic.py --out results/logistic --seed 20261105
python scripts/run_markov_simulation.py --out results/markov --seed 20261105
```

The complete logistic sweep is computationally substantial. For a rapid check, use the short settings documented in [reproducibility](docs/reproducibility.md). Scripts write outputs under `results/` (excluded from Git). Additional figure-generation scripts are available for the manuscript-specific artwork, but they are **not necessary** to use OPS/OFS.

**Scientific provenance:** the Markov certificate in `ordfisher/data/markov_q_exact.json` is a **new, exact rational construction** satisfying the same stated mathematical constraints as the paper's two examples. It is **not a recovered copy** of the historical file. The logistic and simulated Markov results use **new, recorded seeds**. Thus the scripts reproduce the experimental protocol and mathematical claims, not necessarily the same historical simulation realizations or image pixels.

## Documentation

- [Mathematics and conventions](docs/mathematics.md)
- [API reference and edge cases](docs/api.md)
- [Reproducing the paper](docs/reproducibility.md)
- [Quickstart notebook](notebooks/Quickstart_OrdFisher.ipynb)
- [Public release checklist](docs/release-checklist.md)
- [Contributing](CONTRIBUTING.md)
- [Release history](CHANGELOG.md)
- [AI assistance, scope, and validation](docs/ai-use.md)

## License and citation

Released source code is prepared under the **BSD 3-Clause License** (see [`LICENSE`](LICENSE)). Software citation metadata is supplied in [`CITATION.cff`](CITATION.cff). The public **v1.0.0 release is archived on Zenodo**: [DOI: 10.5281/zenodo.24139412](https://doi.org/10.5281/zenodo.24139412). When using the package, cite this versioned software release and the relevant original OPS/OFS scientific papers.

The software authors are Giovanni Campanini, Fernando Pose, Carlos Garcia, and Francisco O. Redelico. The maintainer is **Francisco O. Redelico**. See [`AUTHORS.md`](AUTHORS.md), and the [AI-assistance disclosure](docs/ai-use.md). Please check the [`LICENSE`](LICENSE) and [`AUTHORS.md`](AUTHORS.md) for the published software terms and author attributions.

## Independent ordinal-count validation

An external Colab execution reported **18/18 exact agreements** between `ordfisher` 1.0.0 and `ordpy` 1.2.3 (Python 3.13.18), covering Gaussian noise, Gaussian AR(1), and logistic map for `m=3,4,5` and `tau=1,2`. Maximum reported probability-vector absolute difference was zero in every case; the `m=4` OPS/OFS/Fisher consistency check also passed. This is a validation of ordinal distributions for the tested cases, not a guarantee of correctness for every input. See [the validation record](docs/ordpy-concordance.md).
