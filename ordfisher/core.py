"""Ordinal Power Spectrum (OPS) and Ordinal Fisher Spectrum (OFS).

Conventions used in Campanini et al.'s OFS manuscript (Nov 2026):
* patterns are argsort (time indices ranked by observation value),
* probability vectors are in lexicographic order on S_m,
* Cayley neighbors differ by exchanging consecutive RANKS (positions in
  the permutation tuple), not by swapping neighboring TIMES,
* L = I - A/(m-1) is the normalized Cayley graph Laplacian,
* spectral powers group all eigenvectors with the SAME eigenvalue,
* OPS power is ||Pi_nu(p-u)||^2, OFS power is ||Pi_nu sqrt(p)||^2,
* Fisher energy is sqrt(p).T @ L @ sqrt(p).
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import permutations
from math import factorial

import numpy as np
from scipy.linalg import eigh


@dataclass(frozen=True)
class CayleyGeometry:
    m: int
    patterns: np.ndarray
    edges: np.ndarray
    laplacian: np.ndarray
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray
    mode_groups: tuple[np.ndarray, ...]
    frequencies: np.ndarray


@lru_cache(maxsize=8)
def cayley_geometry(m: int) -> CayleyGeometry:
    """Return the adjacent-RANK Cayley graph, spectral decomposition.

    The full OPS/OFS needs every eigenvalue and eigenspace: dense diagonalization
    is used, with a documented limit m<=6 (N=m! <= 720). No silently truncated
    sparse eigenspectrum is ever substituted for the full spectrum.
    """
    if isinstance(m, bool) or not isinstance(m, (int, np.integer)) or not 2 <= m <= 6:
        raise ValueError("Full Cayley OPS/OFS supports integer 2 <= m <= 6 (dense eigen-decomposition).")
    m = int(m)
    pattern_tuples = tuple(permutations(range(m)))
    pats = np.array(pattern_tuples, dtype=np.int16)
    ix = {p: i for i, p in enumerate(pattern_tuples)}
    n = factorial(m)
    d = m - 1
    adjacency = np.zeros((n, n), dtype=np.float64)
    edges: list[tuple[int, int]] = []
    for i, p in enumerate(pattern_tuples):
        for k in range(m - 1):
            nxt = list(p)
            nxt[k], nxt[k+1] = nxt[k+1], nxt[k]
            j = ix[tuple(nxt)]
            adjacency[i,j] = 1
            if i < j:
                edges.append((i, j))
    assert np.array_equal(adjacency, adjacency.T)
    assert np.all(adjacency.sum(axis=1) == d)
    L = np.eye(n) - adjacency / d
    w, vec = eigh(L, check_finite=False)
    # Numerically equal eigenvalues have basis-dependent eigenvectors; group
    # entire degenerate eigenspaces to make spectral powers basis-independent.
    idx_groups = np.split(np.arange(n), np.flatnonzero(np.diff(w)>1e-8)+1)
    freqs = np.array([np.mean(w[g]) for g in idx_groups])
    assert abs(freqs[0]) < 1e-10
    freqs[0] = 0.0
    assert np.isclose(freqs[-1], 2.0, atol=1e-8)
    edges_array = np.array(edges, dtype=np.int32)
    for array in (pats, edges_array, L, w, vec, freqs, *idx_groups):
        array.setflags(write=False)
    return CayleyGeometry(m, pats, edges_array, L,
                          w, vec, tuple(idx_groups), freqs)


def _validate_probabilities(p, *, tol: float=1e-9) -> np.ndarray:
    p = np.asarray(p, dtype=np.float64)
    if p.ndim != 1:
        raise ValueError("Expected one-dimensional probability vector.")
    if not np.all(np.isfinite(p)) or np.any(p < 0):
        raise ValueError("Probabilities must be finite and nonnegative.")
    if not np.isclose(p.sum(), 1.0, rtol=0, atol=tol):
        raise ValueError(f"Probabilities sum to {p.sum()}, not 1. Normalization is not automatic.")
    return p


def lehmer_indices(patterns: np.ndarray) -> np.ndarray:
    """Map argsort permutation rows to indices in lexicographic S_m order."""
    arr = np.asarray(patterns, dtype=np.int64)
    if arr.ndim != 2 or arr.shape[1] < 2:
        raise ValueError("Expected an array of pattern rows, shape (n_windows,m).")
    m = arr.shape[1]
    # Important for finite-sample speed: all operations are vectorized over windows.
    index = np.zeros(len(arr),dtype=np.int64)
    for i in range(m-1):
        index += np.sum(arr[:,i+1:] < arr[:,[i]],axis=1)*factorial(m-i-1)
    return index


def ordinal_probabilities(data, *, m: int=4, tau: int=1,
                          engine: str="numpy", ties: str="stable",
                          return_diagnostics: bool=False):
    """Calculate every one of the m! probabilities, including missing patterns.

    data: 1D real-valued series. Window length is 1+(m-1)*tau.
    engine: 'numpy' (fast) or 'ordpy' (optional installed dependency).
    ties: 'stable' sorts equal observations by time index; 'raise' rejects them.
    Raw probability vector always follows lexicographic permutation order.
    """
    # Counting ordinal patterns must not trigger a dense Cayley eigendecomposition.
    if isinstance(m, bool) or not isinstance(m, (int, np.integer)) or not 2 <= m <= 6:
        raise ValueError("m must be an integer in [2, 6]")
    m = int(m)
    if not isinstance(tau, (int,np.integer)) or isinstance(tau,bool) or tau < 1:
        raise ValueError("tau must be a positive integer")
    tau=int(tau)
    data=np.asarray(data,dtype=np.float64)
    if data.ndim!=1 or not np.all(np.isfinite(data)):
        raise ValueError("data must be a finite, one-dimensional array")
    windows_n = data.size - (m-1)*tau
    if windows_n < 1:
        raise ValueError("Insufficient observations for m and tau")
    if ties not in ("stable", "raise"):
        raise ValueError("ties must be 'stable' or 'raise'")
    if engine not in ("numpy", "ordpy"):
        raise ValueError("engine must be 'numpy' or 'ordpy'")
    from numpy.lib.stride_tricks import sliding_window_view
    windows = sliding_window_view(data, (m-1)*tau+1)[:,::tau]
    # Counting *windows* with any tie, not total pairs of ties.
    sorted_windows=np.sort(windows,axis=1)
    tie_windows=int(np.count_nonzero(np.any(np.diff(sorted_windows,axis=1)==0,axis=1)))
    if tie_windows and ties=="raise":
        raise ValueError(f"Found {tie_windows} tied ordinal windows (set ties='stable' to keep them)")
    if engine=="numpy":
        order = np.argsort(windows,axis=1,kind="stable")
        counts=np.bincount(lehmer_indices(order),minlength=factorial(m))
        p=counts.astype(float)/windows_n
    else:
        if tie_windows:
            raise ValueError(
                "ordpy tie handling may differ from the stable NumPy policy. "
                "Use engine='numpy' for tied data or resolve ties before using ordpy."
            )
        try:
            import ordpy
        except ImportError as exc:
            raise ImportError("Optional dependency missing: pip install 'ordfisher[ordpy]'") from exc
        # ordpy returns only the supported patterns unless return_missing=True;
        # even when missing patterns are included, ordpy ordering need not be
        # lexicographic. We explicitly reorder by the same Lehmer index.
        pp, probs = ordpy.ordinal_distribution(data, dx=m, taux=tau, return_missing=True)
        pp=np.asarray(pp,dtype=int)
        probs=np.asarray(probs,dtype=float)
        if pp.shape != (factorial(m),m) or probs.shape != (factorial(m),):
            raise RuntimeError("ordpy returned an unexpected ordinal distribution shape")
        p=np.zeros(factorial(m),dtype=float)
        idx=lehmer_indices(pp)
        if len(np.unique(idx))!=factorial(m):
            raise RuntimeError("ordpy returned duplicated or missing permutation labels")
        p[idx]=probs
    diag={"n_observations":len(data), "n_windows":windows_n,
          "n_tied_windows":tie_windows,"m":int(m),"tau":tau,"engine":engine}
    return (p,diag) if return_diagnostics else p


def _spectral_power(geom:CayleyGeometry, signal:np.ndarray) -> np.ndarray:
    coeff=geom.eigenvectors.T@signal
    return np.array([float(np.sum(coeff[g]**2)) for g in geom.mode_groups])


@dataclass(frozen=True)
class Spectrum:
    m: int
    patterns: np.ndarray
    probabilities: np.ndarray
    frequencies: np.ndarray
    ofs: np.ndarray
    ops: np.ndarray
    fisher_energy: float
    fisher_profile: np.ndarray | None
    ofs_profile: np.ndarray | None
    ops_profile: np.ndarray | None
    s0_ofs: float
    s0_ops: float
    shannon_entropy: float
    renyi_half_entropy: float

    def ops_ofs_tv(self) -> float:
        if self.ops_profile is None or self.ofs_profile is None:
            raise ValueError("Normalized profiles undefined at the uniform distribution")
        return total_variation(self.ops_profile,self.ofs_profile)

    def centroids(self) -> tuple[float,float]:
        if self.ops_profile is None or self.ofs_profile is None:
            raise ValueError("Normalized centroids undefined at the uniform distribution")
        nu=self.frequencies[1:]
        return float(nu@self.ops_profile), float(nu@self.ofs_profile)


def total_variation(a, b)->float:
    x=np.asarray(a,dtype=float);y=np.asarray(b,dtype=float)
    if x.shape!=y.shape or x.ndim!=1:
        raise ValueError("TV distance requires same-shape 1D profiles")
    return float(np.sum(np.abs(x-y))/2)


def spectra_from_probabilities(p, *, m: int=4) -> Spectrum:
    """Exact *definition* from the manuscript, numerically evaluated by eigh.

    Positive-frequency normalized OFS/OPS shapes have sum 1. At the uniform
    distribution both profiles and spectral centroids are undefined (None).
    """
    geom=cayley_geometry(m)
    p=_validate_probabilities(p)
    if p.shape != (factorial(m),):
        raise ValueError(f"Expected {factorial(m)} lexicographically ordered probabilities")
    q=np.sqrt(p)
    u=np.full(len(p),1/len(p))
    f=_spectral_power(geom,q)
    power=_spectral_power(geom,p-u)
    # The connected graph has a constant zero eigenvector. Its OPS power
    # is mathematically zero, and its OFS power is the squared affinity to u.
    power[0]=0.0
    f[0]=float((np.sum(q)**2)/len(p))
    # A numerical eigensolver gives eps-sized positive-frequency powers for
    # uniform p; we explicitly detect the mathematical uniform case.
    uniform=np.all(p==u) or np.max(np.abs(p-u)) < 5e-14
    if uniform:
        f[1:]=0.0
        power[:]=0.0
    s0=float(f[1:].sum());sP=float(power[1:].sum())
    energy=float(np.dot(geom.frequencies[1:],f[1:]))
    if abs(energy)<1e-14: energy=0.0
    ofs_profile=None if s0<1e-14 else f[1:]/s0
    ops_profile=None if sP<1e-14 else power[1:]/sP
    fisher_profile=None if energy<1e-14 else geom.frequencies[1:]*f[1:]/energy
    nonzero=p[p>0]
    entropy=-float(np.sum(nonzero*np.log(nonzero)))/np.log(len(p))
    renyi=float(2*np.log(np.sum(q))/np.log(len(p)))
    return Spectrum(m,geom.patterns,p,geom.frequencies,f,power,energy,
                    fisher_profile,ofs_profile,ops_profile,s0,sP,entropy,renyi)


def spectra_from_series(data, *, m:int=4,tau:int=1,engine:str="numpy",ties:str="stable") -> Spectrum:
    return spectra_from_probabilities(ordinal_probabilities(data,m=m,tau=tau,engine=engine,ties=ties),m=m)


def check_spectral_identities(p, *, m:int=4, atol:float=5e-10)->dict[str,float]:
    """Recompute required identities directly from p and graph edges."""
    geom=cayley_geometry(m)
    r=spectra_from_probabilities(p,m=m)
    p=r.probabilities
    q=np.sqrt(p)
    n=len(p)
    edges=geom.edges
    energy_edges=float(np.sum((q[edges[:,0]]-q[edges[:,1]])**2)/(m-1))
    s0_renyi=1-(np.sum(q)**2)/n
    s0_ops_euclidean=float(np.sum((p-1/n)**2))
    d={
        "fisher_moment_minus_edges":r.fisher_energy-energy_edges,
        "fisher_moment_minus_quadratic_form":r.fisher_energy-float(q@geom.laplacian@q),
        "ofs_sum_minus_one":float(r.ofs.sum())-1.0,
        "ofs_positive_minus_renyi":r.s0_ofs-s0_renyi,
        "ops_positive_minus_euclidean":r.s0_ops-s0_ops_euclidean,
    }
    assert all(abs(v)<atol for v in d.values()), d
    return d


def keller_order_m3() -> tuple[tuple[int,...],...]:
    return ((0,1,2),(0,2,1),(2,0,1),(1,0,2),(1,2,0),(2,1,0))


def path_energy(p, order, *, m:int)->float:
    """E_rho=1/2*sum (sqrt(p_{rho[k+1]})-sqrt(p_{rho[k]}))^2."""
    p=_validate_probabilities(p)
    if p.shape != (factorial(m),):raise ValueError("Invalid p length")
    ix={tuple(pattern): i for i,pattern in enumerate(cayley_geometry(m).patterns)}
    ids=[ix[tuple(v)] for v in order]
    if len(set(ids))!=factorial(m):raise ValueError("order must list all patterns exactly once")
    return float(0.5*np.sum(np.diff(np.sqrt(p)[ids])**2))
