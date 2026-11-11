#!/usr/bin/env python3
"""Minimal example for arbitrary time series (without ordpy dependency)."""
import numpy as np
from ordfisher import spectra_from_series, check_spectral_identities

rng=np.random.default_rng(42)
series=rng.normal(size=10000)
result=spectra_from_series(series,m=4,tau=1)
print('Graph frequencies:',result.frequencies)
print('OPS powers:',result.ops)
print('OFS powers:',result.ofs)
print('Scalar Fisher energy:',result.fisher_energy)
print('OPS/OFS normalized distance:',result.ops_ofs_tv())
print('Identity residuals:',check_spectral_identities(result.probabilities,m=4))
