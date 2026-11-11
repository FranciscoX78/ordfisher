# Mathematical definitions and conventions

## Ordinal patterns

A time series `x` is represented by length-`m` windows sampled at embedding delay `tau`. Each window maps to the permutation of time indices `0,...,m-1` that sorts the observations from smallest to largest. In other words, the pattern is the *argsort* permutation, not a vector of ranks attached to the original time positions.

The probability vector `p` has `N=m!` entries in **lexicographic** permutation order, including entries equal to zero. All subsequent matrices use that identical ordering.

When values tie exactly, the NumPy backend uses a stable sort (earlier time index first), unless `ties='raise'` is requested. The optional `ordpy` backend is reserved for tie-free input.

## Adjacent-rank graph

The undirected graph has one vertex per permutation. Its edges connect patterns differing by a swap of **adjacent positions in the argsort permutation** (adjacent *ranks*), not a swap of neighboring time indices. Let `A` be the adjacency matrix and let `d=m-1` be the vertex degree. The normalized Laplacian is

\[
L=I-\frac{A}{d}.
\]

For each **distinct** eigenvalue `nu` of `L`, the orthogonal projector `Pi_nu` denotes projection onto its entire eigenspace. Degenerate eigenvalues are grouped to make OPS and OFS independent of arbitrary choices of eigenvectors within that eigenspace.

## Spectral definitions

Define `u=(1/N,...,1/N)` and `q=sqrt(p)` componentwise. Then

\[
\mathrm{OPS}(\nu;p)=\|\Pi_\nu(p-u)\|^2,\qquad
\mathrm{OFS}(\nu;p)=\|\Pi_\nu q\|^2.
\]

The zero-frequency OFS power is the squared affinity with the uniform distribution; the zero-frequency OPS power is zero by construction.

The Cayley square-root Fisher energy is

\[
I_G(p)=q^\top Lq
 =\frac1{m-1}\sum_{\{i,j\}\in E}(\sqrt{p_i}-\sqrt{p_j})^2
 =\sum_{\nu>0}\nu\,\mathrm{OFS}(\nu;p).
\]

Its positive-frequency allocation is

\[
W_G(\nu;p)=\frac{\nu\,\mathrm{OFS}(\nu;p)}{I_G(p)},\qquad\nu>0.
\]

By contrast, the normalized OPS/OFS profiles divide *unweighted* positive-frequency powers by their respective positive-frequency totals. The distinction is important when interpreting distances between spectra.

For `S0_OFS=sum_{nu>0} OFS(nu)` and `S0_OPS=sum_{nu>0} OPS(nu)`, the identities are

\[
\sum_{\nu}\mathrm{OFS}(\nu;p)=1,
\quad S0_{\rm OFS}=1-\frac{(\sum_i\sqrt{p_i})^2}{N},
\quad S0_{\rm OPS}=\sum_i(p_i-1/N)^2.
\]

`S0_OFS` is exactly determined by Rényi permutation entropy of order one-half. **It is not a new independent scalar**. The nontrivial information is in the *resolved* spectrum and its Fisher-energy weighting. The mean of path-based Fisher energies over all `(m!)!` numerical enumerations equals `S0_OFS`; that averaging is a global identity, **not** the algorithm used to calculate OFS.

## What the implementation is and is not

- It uses a dense, real symmetric eigendecomposition for the **entire** graph spectrum. `m<=6` is an explicit computational limit, not a theorem about admissible pattern lengths.
- Zero modes and numerically vanishing positive-frequency totals are treated separately; normalized profiles are undefined when the denominator vanishes.
- `fisher_energy` is a **graph gradient energy**, not a parametric Fisher information matrix.
- Graph frequencies describe smooth and rapidly changing structures *on the ordinal graph*, not time-series Fourier frequencies.
- The output of `ordfisher` is a property of the specified marginal ordinal distribution and graph. It does **not**, on its own, reconstruct successive ordinal transitions.
