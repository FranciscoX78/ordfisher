"""Build and visualize the adjacent-rank Cayley graph used by OPS/OFS.

For a length-m ordinal pattern (the indices ordered by increasing value),
two patterns are neighbors precisely when one is obtained from the other
by swapping positions k and k+1 in the permutation (adjacent *ranks*).

The full graph has m! nodes, degree m-1 and m!(m-1)/2 undirected edges.
This module constructs the graph without diagonalizing its Laplacian.
"""
from __future__ import annotations

from itertools import permutations
from math import factorial
from pathlib import Path
from importlib.resources import files
import json
from typing import Literal

import networkx as nx


def cayley_graph(m: int = 4) -> nx.Graph:
    """Return the full adjacent-rank Cayley graph of S_m.

    Node labels are permutations, e.g. (0, 1, 2, 3). An edge joins
    two nodes if and only if their adjacent *ranks* were exchanged.
    Its ``generator`` attribute (1..m-1) denotes the pair of ranks
    exchanged (1-based); the index positions in the tuple are 0-based.

    Supported m: 2..7. Plotting graphs for m>=5 may be crowded.
    This combinatorial function is independent of the m<=6 dense
    eigendecomposition bound of ``cayley_geometry``.
    """
    if isinstance(m, bool) or not isinstance(m, int) or not 2 <= m <= 7:
        raise ValueError("m must be an integer between 2 and 7")

    graph = nx.Graph(m=m, convention="swap adjacent ranks, not adjacent times")
    for p in permutations(range(m)):
        graph.add_node(p)
    for p in graph.nodes:
        for k in range(m - 1):
            q = list(p)
            q[k], q[k + 1] = q[k + 1], q[k]
            q = tuple(q)
            # Adding each undirected edge only once also preserves its generator.
            if p < q:
                graph.add_edge(p, q, generator=k + 1)

    expected_nodes = factorial(m)
    expected_edges = expected_nodes * (m-1) // 2
    assert graph.number_of_nodes() == expected_nodes
    assert graph.number_of_edges() == expected_edges
    return graph


def plot_cayley_graph(
    m: int = 4,
    *,
    layout: Literal["appendix", "spring", "kamada_kawai", "circular"] | None = None,
    color_by_generator: bool = True,
    with_labels: bool | None = None,
    seed: int = 17,
    figsize: tuple[float, float] | None = None,
    save: str | Path | None = None,
    ax=None,
):
    """Draw the full Cayley graph and optionally save it as PNG/PDF/SVG.

    Returns ``(fig, ax)``. The legend names the rank swaps (1,2),
    (2,3), etc.  In ``spring`` layout the ``seed`` fixes positions.
    For m=4, the default layout is "appendix": the *exact original*
    coordinates and visual conventions of Fig. D1 in the Chaos manuscript,
    with the frozen original placement. Use layout="spring" to obtain the
    former automatic rendering. For other m the default remains "spring".
    The appendix preset supports m=4 only. Custom plot options intentionally
    change the appearance from the publication figure.

    For large m, labels are disabled by default to avoid overcrowding.
    The full graph (not an eigenvector projection) is being plotted.
    """
    # Import matplotlib only for users who ask for a drawing.
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    graph = cayley_graph(m)
    if layout is None:
        layout = "appendix" if m == 4 else "spring"
    if with_labels is None:
        with_labels = m <= 4

    if layout == "appendix":
        if m != 4:
            raise ValueError("layout='appendix' is the published S4 figure and requires m=4")
        # This path deliberately uses the original graph insertion order,
        # precisely the same node and edge draw order as the manuscript source.
        # A frozen placement avoids dependency on NetworkX layout algorithms
        # or their version-specific default traversal order.
        original = nx.Graph()
        for p in permutations(range(4)):
            original.add_node(p)
            for j in range(3):
                q = list(p)
                q[j], q[j + 1] = q[j + 1], q[j]
                original.add_edge(p, tuple(q), kind=j + 1)
        assert {frozenset((u,v)) for u,v in original.edges()} == {frozenset((u,v)) for u,v in graph.edges()}
        placements = json.loads(files("ordfisher").joinpath("data", "cayley_s4_appendix_layout.json").read_text())
        pos = {p: tuple(placements["".join(map(str,p))]) for p in original}
        graph = original
        created_ax = ax is None
        if created_ax:
            fig, ax = plt.subplots(figsize=figsize if figsize is not None else (10.2, 8.0))
        else:
            fig = ax.figure
        colors = {1:'#8296AF', 2:'#237C86', 3:'#D28B50'}
        for j in (1,2,3):
            edgelist = [(u,v) for u,v,d in graph.edges(data=True) if d['kind']==j]
            nx.draw_networkx_edges(graph, pos, edgelist=edgelist,
                                   edge_color=colors[j] if color_by_generator else '0.6',
                                   width=1.8, ax=ax)
        # Parity coloring is part of the original appendix figure.
        node_colors = ['#EEF3F8' if sum(p[i] > p[j] for i in range(4) for j in range(i+1,4))%2==0 else '#E5F3F0' for p in graph]
        nx.draw_networkx_nodes(graph, pos, node_size=750, node_color=node_colors,
                               edgecolors='#344D61', linewidths=1.3, ax=ax)
        if with_labels:
            nx.draw_networkx_labels(graph, pos,
                                    labels={p:''.join(map(str,p)) for p in graph},
                                    font_size=7.8, font_weight='medium', ax=ax)
        if color_by_generator:
            ax.legend(handles=[Line2D([0],[0], color=colors[j], lw=2.8,
                    label=f'Swap ranks {j}–{j+1}') for j in (1,2,3)],
                    frameon=False, loc='lower right', fontsize=10)
        ax.set_aspect('equal')
        ax.axis('off')
        if created_ax and figsize is None:
            # Use exactly the same panel margins as the original manuscript script.
            fig.subplots_adjust(left=.04, right=.96, top=.98, bottom=.04)
        if save is not None:
            filename = Path(save)
            if filename.suffix.lower() not in {'.png','.pdf','.svg'}:
                raise ValueError('save must end in .png, .pdf, or .svg')
            filename.parent.mkdir(parents=True, exist_ok=True)
            opts = {'bbox_inches':'tight', 'pad_inches':.08}
            if filename.suffix.lower() == '.png':
                opts['dpi'] = 200
            fig.savefig(filename, **opts)
        return fig, ax

    if layout == "spring":
        pos = nx.spring_layout(graph, seed=seed, iterations=400)
    elif layout == "kamada_kawai":
        pos = nx.kamada_kawai_layout(graph)
    elif layout == "circular":
        pos = nx.circular_layout(graph)
    else:
        raise ValueError("layout must be 'appendix', 'spring', 'kamada_kawai' or 'circular'")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize if figsize is not None else (10,8))
    else:
        fig = ax.figure

    if color_by_generator:
        palette = plt.get_cmap("tab10")
        legend_handles = []
        for generator in range(1, m):
            edges = [(u,v) for u,v,d in graph.edges(data=True)
                     if d["generator"] == generator]
            color = palette((generator - 1) % 10)
            nx.draw_networkx_edges(graph, pos, ax=ax, edgelist=edges,
                                   edge_color=[color], width=1.8, alpha=0.75)
            legend_handles.append(Line2D([0],[0], color=color, linewidth=2,
                                         label=f"Ranks {generator} and {generator+1}"))
    else:
        nx.draw_networkx_edges(graph, pos, ax=ax, edge_color="0.60",
                               width=1.6, alpha=0.85)

    nx.draw_networkx_nodes(graph, pos, ax=ax, node_color="white",
                           edgecolors="#28374A", linewidths=1.4,
                           node_size=440 if m <= 3 else 390 if m == 4 else 80)

    if with_labels:
        if m <= 10:
            labels = {p: "".join(map(str,p)) for p in graph}
        else:
            labels = {p: str(p) for p in graph}
        nx.draw_networkx_labels(graph, pos, labels=labels, ax=ax,
                                font_size=9 if m <= 3 else 7 if m == 4 else 5,
                                font_family="DejaVu Sans")

    ax.set_title(f"Adjacent-rank Cayley graph of $S_{m}$: "
                 f"{graph.number_of_nodes()} patterns, {graph.number_of_edges()} edges",
                 fontsize=12, pad=15)
    ax.set_aspect("equal", adjustable="datalim")
    ax.axis("off")
    if color_by_generator:
        ax.legend(handles=legend_handles, title="Adjacent-rank swaps",
                  loc="upper left", bbox_to_anchor=(1.01,1), frameon=False,
                  fontsize=8, title_fontsize=8)
    fig.tight_layout()
    if save is not None:
        filename = Path(save)
        if filename.suffix.lower() not in {".png",".pdf",".svg"}:
            raise ValueError("save must end in .png, .pdf, or .svg")
        filename.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(filename, dpi=250, bbox_inches="tight")
    return fig, ax
