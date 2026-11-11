#!/usr/bin/env python3
"""Optional cross-check with ordpy on data without ties."""
import importlib.util
if importlib.util.find_spec('ordpy') is None:
    print("SKIPPED: ordpy is optional and not installed. Install: pip install 'ordfisher[ordpy]'")
    raise SystemExit(0)
import numpy as np
from ordfisher import ordinal_probabilities
rng=np.random.default_rng(20261105)
x=rng.standard_normal(2000)
for m in (3,4,5):
    a=ordinal_probabilities(x,m=m,engine='numpy',ties='raise')
    b=ordinal_probabilities(x,m=m,engine='ordpy',ties='raise')
    err=float(np.max(np.abs(a-b)))
    print(f'm={m}: max difference against ordpy = {err:.3g}')
    assert err<1e-12
print('ordpy/Numpy ordinal distributions agree for all tested m.')
