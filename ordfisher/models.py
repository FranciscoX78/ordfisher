"""Synthetic logistic and exact Gaussian AR(1) ordinal distributions."""
from __future__ import annotations
from math import pi
import numpy as np
from scipy.optimize import brentq
from .core import cayley_geometry, spectra_from_probabilities, total_variation


def _logistic_pure(r:float,x0:float,burn:int,retained:int)->np.ndarray:
    x=float(x0)
    for _ in range(burn):
        x=r*x*(1-x)
    out=np.empty(retained,dtype=np.float64)
    for t in range(retained):
        out[t]=x
        x=r*x*(1-x)
    return out


try:
    from numba import njit
    _logistic_fast=njit(cache=True)(_logistic_pure)
    _HAS_NUMBA=True
except ImportError:
    _HAS_NUMBA=False


def logistic_series(r:float,x0:float, *, burn:int=5000,retained:int=120000,
                    use_numba:bool=True)->np.ndarray:
    """Discard `burn` iterates, retain `retained`, from x_{t+1}=r*x_t*(1-x_t)."""
    if not (0<r<=4) or not (0<x0<1):raise ValueError("Need 0<r<=4 and 0<x0<1")
    if burn<0 or retained<1:raise ValueError("Invalid burn/retained lengths")
    fn=_logistic_fast if use_numba and _HAS_NUMBA else _logistic_pure
    return fn(float(r),float(x0),int(burn),int(retained))


def ar1_gaussian_ordinal_probabilities(phi:float,*,m:int=4)->np.ndarray:
    """Exact (up to floating point evaluation) Gaussian AR(1) ordinal p for m=4.

    Formula: orthant probability of three joint Gaussian differences,
    1/8 + [asin(R12)+asin(R13)+asin(R23)]/(4*pi).
    Covariance matrix Cij=phi**abs(i-j).
    """
    if m!=4:raise ValueError("The analytic arcsine formula implemented here is for m=4 only")
    if not -1<phi<1:raise ValueError("Stationary Gaussian AR(1) requires |phi|<1")
    geom=cayley_geometry(4)
    row=np.arange(4)
    cov=np.power(float(phi),np.abs(row[:,None]-row[None,:]))
    p=np.empty(24)
    for t,perm in enumerate(geom.patterns):
        d=np.zeros((3,4))
        for j in range(3):
            d[j,perm[j+1]]=1
            d[j,perm[j]]=-1
        V=d@cov@d.T
        denom=np.sqrt(np.outer(V.diagonal(),V.diagonal()))
        R=np.clip(V/denom,-1,1)
        p[t]=0.125+sum(np.arcsin(R[i,j]) for i,j in [(0,1),(0,2),(1,2)])/(4*pi)
    assert np.isclose(np.sum(p),1,atol=1e-10)
    if np.min(p)<-1e-12:raise RuntimeError("Unphysical probability in AR(1) computation")
    p=np.maximum(p,0)
    return p/p.sum()


def equal_fisher_ar1_pair()->dict:
    """Find phi*>0 having IG(phi*)=IG(-0.5), then compare energy profiles."""
    target=spectra_from_probabilities(ar1_gaussian_ordinal_probabilities(-0.5),m=4)
    def resid(phi):
        return spectra_from_probabilities(ar1_gaussian_ordinal_probabilities(phi),m=4).fisher_energy-target.fisher_energy
    star=float(brentq(resid,.65,.66,xtol=1e-14))
    other=spectra_from_probabilities(ar1_gaussian_ordinal_probabilities(star),m=4)
    return {"phi_negative":-0.5,"phi_positive":star,
            "fisher_energy_negative":target.fisher_energy,
            "fisher_energy_positive":other.fisher_energy,
            "tv_fisher_energy_profiles":total_variation(target.fisher_profile,other.fisher_profile),
            "tv_normalized_ofs_profiles":total_variation(target.ofs_profile,other.ofs_profile)}
