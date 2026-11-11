"""Public API and packaging-oriented regression tests."""
import importlib.resources as resources
import json
import numpy as np
import pytest
from ordfisher import (
    __version__, cayley_geometry, ordinal_probabilities,
    spectra_from_probabilities, cayley_graph,
)


def test_release_version():
    assert __version__ == '1.0.0'


def test_graph_cache_protected():
    geom = cayley_geometry(4)
    with pytest.raises(ValueError):
        geom.laplacian[0,0] = 42
    with pytest.raises(ValueError):
        geom.patterns[0,0] = 3
    with pytest.raises(ValueError):
        geom.frequencies[0] = 4


def test_input_validation():
    for p in ([1/6]*5, [-0.2]+[0.24]*5, [float('nan')]+[0.2]*5):
        with pytest.raises(ValueError):
            spectra_from_probabilities(p, m=3)
    for m in (0, 7, 2.5, True):
        with pytest.raises(ValueError):
            ordinal_probabilities([1,2,3,4,5],m=m)


def test_optional_ordpy_adapter_reorders_and_handles_ties(monkeypatch):
    # Fake ordpy verifies our adapter's relabeling and fail-safe tie policy.
    import sys, types
    expected = ordinal_probabilities([5,3,2,7,1,8,4],m=3)
    perms = np.asarray(list(__import__('itertools').permutations(range(3))))
    permuted = np.asarray([5,3,0,1,4,2])
    fake = types.SimpleNamespace(ordinal_distribution=lambda *args, **kwargs: (perms[permuted], expected[permuted]))
    monkeypatch.setitem(sys.modules, 'ordpy', fake)
    actual = ordinal_probabilities([5,3,2,7,1,8,4],m=3,engine='ordpy')
    assert np.allclose(actual,expected)
    with pytest.raises(ValueError, match='tie handling'):
        ordinal_probabilities([0,0,1,2,3],m=3,engine='ordpy',ties='stable')


def test_data_resources_in_package():
    package_data = resources.files('ordfisher').joinpath('data')
    obj = json.loads(package_data.joinpath('markov_q_exact.json').read_text())
    layout = json.loads(package_data.joinpath('cayley_s4_appendix_layout.json').read_text())
    assert obj['bins'] == 8
    assert len(obj['plus']['numerators']) == 4096
    assert len(layout) == 24


def test_graph_degree_and_modes_m2():
    g=cayley_graph(2)
    assert (g.number_of_nodes(),g.number_of_edges()) == (2,1)
    spec=spectra_from_probabilities([0.7,0.3],m=2)
    assert spec.frequencies.shape == (2,)
    assert spec.ofs_profile.shape == (1,)
