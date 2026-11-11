import numpy as np
import pytest
from math import factorial
from ordfisher import (cayley_geometry,ordinal_probabilities,spectra_from_probabilities,
                       spectra_from_series,check_spectral_identities, keller_order_m3,path_energy)
from ordfisher.core import lehmer_indices

@pytest.mark.parametrize('m',[2,3,4,5])
def test_graph_structure(m):
    g=cayley_geometry(m)
    assert len(g.patterns)==factorial(m)
    assert len(g.edges)==factorial(m)*(m-1)//2
    assert abs(g.frequencies[0])<1e-12
    assert abs(g.frequencies[-1]-2)<1e-10
    assert np.allclose(g.laplacian,g.laplacian.T)

@pytest.mark.parametrize('m',[3,4,5])
def test_identities_random(m):
    rng=np.random.default_rng(m)
    p=rng.dirichlet(np.ones(factorial(m)))
    check_spectral_identities(p,m=m)
    r=spectra_from_probabilities(p,m=m)
    assert abs(r.ofs.sum()-1)<1e-12
    assert abs(sum(r.ops_profile)-1)<1e-12
    assert abs(sum(r.ofs_profile)-1)<1e-12
    assert abs(sum(r.fisher_profile)-1)<1e-12


def test_uniform_is_undefined_for_normalized_profiles():
    r=spectra_from_probabilities(np.ones(24)/24,m=4)
    assert r.ofs_profile is None and r.ops_profile is None
    assert r.fisher_profile is None
    assert r.fisher_energy==0
    with pytest.raises(ValueError):r.ops_ofs_tv()


def test_missing_patterns_are_not_dropped():
    p=np.zeros(24);p[0]=1
    r=spectra_from_probabilities(p,m=4)
    check_spectral_identities(p,m=4)
    assert r.s0_ofs>0


def test_numpy_bandt_pompe():
    p=ordinal_probabilities([4,7,9,10,6,11,3],m=2)
    assert np.allclose(p,[2/3,1/3])
    rng=np.random.default_rng(102)
    x=rng.standard_normal(1000)
    for m in (3,4):
        counts=ordinal_probabilities(x,m=m)
        order=np.argsort(np.lib.stride_tricks.sliding_window_view(x,m),axis=1,kind='stable')
        from collections import Counter
        counter=Counter(map(tuple,order))
        perms=list(map(tuple,cayley_geometry(m).patterns))
        assert np.allclose(counts,np.array([counter.get(perm,0) for perm in perms])/(len(x)-m+1))


def test_ties_are_explicit():
    x=[0,0,1,1,2,2]
    with pytest.raises(ValueError,match='tied'): ordinal_probabilities(x,m=3,ties='raise')
    p,diag=ordinal_probabilities(x,m=3,ties='stable',return_diagnostics=True)
    assert diag['n_tied_windows']>0
    assert abs(p.sum()-1)<1e-12


def test_keller_path_energy():
    p=np.array([1/3,1/15,2/15,1/5,4/15,0])
    assert abs(path_energy(p,keller_order_m3(),m=3)-0.199067676998)<1e-10


def test_input_errors():
    with pytest.raises(ValueError):spectra_from_probabilities([1/2,1/2],m=3)
    with pytest.raises(ValueError):cayley_geometry(7)
    with pytest.raises(ValueError):ordinal_probabilities([1,2],m=4)
    with pytest.raises(ValueError):ordinal_probabilities([1,np.nan,2,3],m=3)
