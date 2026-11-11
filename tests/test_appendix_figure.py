"""Regression tests for the exact original figure used in the Chaos appendix."""
from pathlib import Path
import os
import json
from importlib.resources import files

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image, ImageChops
import networkx as nx
import pytest

from ordfisher import cayley_graph, plot_cayley_graph


def test_appendix_fixed_coordinates_and_graph():
    d = json.loads(files('ordfisher').joinpath('data', 'cayley_s4_appendix_layout.json').read_text())
    assert len(d) == 24
    assert set(d) == {''.join(map(str,p)) for p in cayley_graph(4)}
    assert len(set(tuple(x) for x in d.values())) == 24
    assert all(len(x) == 2 for x in d.values())
    g = cayley_graph(4)
    assert g.number_of_edges() == 36
    assert nx.is_planar(g)


def test_appendix_default_and_explicit_match(tmp_path):
    a = tmp_path/'default.png'
    b = tmp_path/'explicit.png'
    fa, _ = plot_cayley_graph(4,save=a)
    fb, _ = plot_cayley_graph(4,layout='appendix',save=b)
    plt.close(fa);plt.close(fb)
    assert a.read_bytes() == b.read_bytes()
    assert Image.open(a).width > 1000
    assert Image.open(a).height > 1000


def test_only_s4_can_use_appendix():
    with pytest.raises(ValueError, match='requires m=4'):
        plot_cayley_graph(3,layout='appendix')


def test_original_pixel_identity_opt_in(tmp_path):
    # Strict identical pixels need the same fonts and Matplotlib rasterizer.
    if os.getenv('ORDFISHER_STRICT_PLOT_PIXELS') != '1':
        pytest.skip('Set ORDFISHER_STRICT_PLOT_PIXELS=1 to compare with the original PNG')
    new = tmp_path/'recreated.png'
    fig, _ = plot_cayley_graph(4,layout='appendix',save=new)
    plt.close(fig)
    reference = Path(__file__).parent/'fixtures'/'figD1_Cayley_m4_original.png'
    aa=Image.open(reference).convert('RGBA')
    bb=Image.open(new).convert('RGBA')
    assert aa.size == bb.size
    assert ImageChops.difference(aa,bb).getbbox() is None
