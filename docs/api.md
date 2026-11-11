# Public API

All functions documented here are importable from `ordfisher` unless otherwise indicated.

## Main functions

| Function | Purpose |
|---|---|
| `ordinal_probabilities(data, m=4, tau=1, engine='numpy', ties='stable')` | Returns a complete lexicographically ordered length-`m!` probability vector, including zeros. |
| `spectra_from_probabilities(p, m=4)` | Computes complete Cayley OPS/OFS spectra from an already-normalized ordinal distribution. |
| `spectra_from_series(data, m=4, tau=1, engine='numpy', ties='stable')` | Computes the ordinal distribution and both spectra from a real-valued time series. |
| `cayley_geometry(m)` | Cached complete graph, edges, Laplacian, eigendecomposition, and distinct-frequency groups (supports `2<=m<=6`). |
| `cayley_graph(m=4)` | NetworkX graph with permutation-valued vertices and edges for adjacent-rank swaps (supports `2<=m<=7`). |
| `plot_cayley_graph(m=4, layout=None, save=None, ...)` | Draws a Cayley graph; returns Matplotlib `(figure, axes)`. |
| `check_spectral_identities(p, m=4)` | Independently compares the Fisher, OFS, OPS and Rényi-related sum identities. |
| `total_variation(a,b)` | Half the sum of the absolute entrywise differences of equal-length vectors. |

## `Spectrum` result fields

A `Spectrum` returned by `spectra_from_series` or `spectra_from_probabilities` exposes:

- `m`, `patterns`, `probabilities`: embedding length, lexicographic permutation labels, and ordinal distribution.
- `frequencies`: **distinct** graph eigenvalues, including zero.
- `ops`, `ofs`: complete spectral power arrays, including the zero mode.
- `fisher_energy`: scalar Cayley Fisher energy.
- `ops_profile`, `ofs_profile`: powers at **positive** frequencies normalized to sum one; `None` at uniformity.
- `fisher_profile`: shares of *frequency-weighted* OFS energy; `None` when Fisher energy vanishes.
- `s0_ops`, `s0_ofs`: positive-frequency unnormalized sums.
- `shannon_entropy`, `renyi_half_entropy`: normalized ordinal entropy values.
- `ops_ofs_tv()`: total-variation distance between the normalized OFS/OPS profiles.
- `centroids()`: the graph-frequency centroids of the normalized OPS and OFS profiles (in that order).

## Input requirements and boundaries

- `m` is an integer between `2` and `6` for complete spectra. Full eigendecomposition costs grow rapidly as `m!` increases.
- `tau` is a strictly positive integer. The input data must be finite, one-dimensional, and long enough for a window of size `1+(m-1)*tau`.
- `p` must be a one-dimensional vector of length exactly `m!`, finite and nonnegative, and sum to one within a strict tolerance. It is **not automatically renormalized**.
- Lexicographic `argsort` ordering is essential: use an explicit conversion if another package supplies a different pattern order.
- Exact ties are rejected by `ties='raise'` or resolved stably by `ties='stable'` for the NumPy engine. When using `engine='ordpy'`, **tie-free input is required**.
- `plot_cayley_graph(4)` uses the manuscript's fixed `appendix` layout. Other layouts are `spring`, `kamada_kawai`, and `circular`. Save formats are PNG, PDF, and SVG.
- **Do not** interpret normalized profiles as zero vectors at uniformity: their denominators vanish. The library returns `None` and methods requiring normalized profiles raise `ValueError`.

## Reuse with ordpy

After installing the optional dependency, a tie-free comparison is:

```python
from ordfisher import ordinal_probabilities

p_a = ordinal_probabilities(x, m=4, engine='numpy', ties='raise')
p_b = ordinal_probabilities(x, m=4, engine='ordpy', ties='raise')
```

An optional independent check is supplied by `scripts/check_ordpy_agreement.py`. If `ordpy` is not installed, the script currently skips its check and prints a message. Do not present a skipped check as a passed independent validation.
