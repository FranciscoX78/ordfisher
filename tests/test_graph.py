"""Tests of the newly exposed Cayley graph and plotting API."""
from math import factorial
from itertools import combinations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import pytest

from ordfisher import cayley_geometry, cayley_graph, plot_cayley_graph


def kendall(p, q):
    positions = {v: i for i, v in enumerate(q)}
    return sum(positions[p[i]] > positions[p[j]]
               for i, j in combinations(range(len(p)), 2))


@pytest.mark.parametrize("m", [2, 3, 4, 5])
def test_graph_counts_degrees_and_adjacency(m):
    graph = cayley_graph(m)
    assert len(graph) == factorial(m)
    assert graph.number_of_edges() == factorial(m)*(m-1)//2
    assert nx.is_connected(graph)
    assert all(d == m-1 for _, d in graph.degree())
    assert all(kendall(p,q) == 1 for p,q in graph.edges())
    assert all(1 <= data["generator"] <= m-1 for _,_,data in graph.edges(data=True))


def test_s3_is_six_cycle():
    assert nx.is_isomorphic(cayley_graph(3), nx.cycle_graph(6))


def test_s4_matches_spectral_geometry():
    g = cayley_graph(4)
    geom = cayley_geometry(4)
    graph_edges = {frozenset((u,v)) for u,v in g.edges()}
    expected_edges = {frozenset((tuple(geom.patterns[a]), tuple(geom.patterns[b])))
                      for a,b in geom.edges}
    assert graph_edges == expected_edges


def test_distance_keller_example():
    g = cayley_graph(3)
    assert nx.shortest_path_length(g,(2,0,1),(1,0,2)) == 3


@pytest.mark.parametrize("m", [3,4])
def test_draw_and_export(m,tmp_path):
    for ext in ("png", "pdf"):
        filename = tmp_path/f"cayley_s{m}.{ext}"
        fig, ax = plot_cayley_graph(m,save=filename)
        assert filename.exists() and filename.stat().st_size > 1000
        assert ax.figure is fig
        plt.close(fig)


def test_input_errors(tmp_path):
    for value in (1, 8, 2.5, True):
        with pytest.raises(ValueError):
            cayley_graph(value)
    with pytest.raises(ValueError,match="layout"):
        plot_cayley_graph(3,layout="bad")
    with pytest.raises(ValueError,match="save"):
        plot_cayley_graph(3,save=tmp_path/'image.bmp')
