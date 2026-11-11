import numpy as np
from ordfisher.markov import (verify_exact_certificate,exact_Q,theoretical_probabilities,
                               ordinal_bin_coefficient_matrix,ordinal_transition_probability,
                               simulate_markov)
from ordfisher import spectra_from_probabilities,total_variation
from ordfisher.models import equal_fisher_ar1_pair,ar1_gaussian_ordinal_probabilities


def test_rational_certificates():
    info=verify_exact_certificate()
    assert info['plus']['max_exact_ordinal_residual']==0
    assert info['minus']['max_exact_ordinal_residual']==0
    for name in ('plus','minus'):
        num,den=exact_Q(name)
        assert np.min(num)>0 and np.sum(num)==den
        assert np.all(num.sum(axis=3)==den//512)
        assert np.all(num.sum(axis=0)==den//512)


def test_exact_transition_probabilities():
    assert str(ordinal_transition_probability(1,'021','120'))=='56/125'
    assert str(ordinal_transition_probability(-1,'021','120'))=='13/250'


def test_markov_ofs_separation():
    a=spectra_from_probabilities(theoretical_probabilities(1),m=4)
    b=spectra_from_probabilities(theoretical_probabilities(-1),m=4)
    assert np.allclose(a.ops,b.ops,atol=1e-12)
    assert abs(a.fisher_energy-b.fisher_energy)<1e-12
    assert abs(total_variation(a.ofs_profile,b.ofs_profile)-0.2312292799484)<1e-10


def test_simulation_is_stationary_in_sample():
    x=simulate_markov('plus',10000,seed=111)
    assert x.shape==(10000,)
    assert np.all((x>=0)&(x<1))
    assert abs(x.mean()-0.5)<.03


def test_ar1_exact_pair():
    z=equal_fisher_ar1_pair()
    assert abs(z['phi_positive']-.6536965296)<1e-8
    assert abs(z['tv_fisher_energy_profiles']-.412758422)<1e-8
    assert abs(z['tv_normalized_ofs_profiles']-.322300900)<1e-8
    assert np.allclose(ar1_gaussian_ordinal_probabilities(0),1/24,atol=1e-14)
