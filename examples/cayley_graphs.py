"""Generate the S4 Cayley figure with the paper's fixed appendix layout."""
import matplotlib.pyplot as plt
from ordfisher import cayley_graph, plot_cayley_graph

g = cayley_graph(4)
assert g.number_of_nodes() == 24
assert g.number_of_edges() == 36
fig, ax = plot_cayley_graph(4, layout='appendix', save='cayley_S4.pdf')
plt.close(fig)
