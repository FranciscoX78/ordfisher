"""Cayley-graph OPS/OFS computations and reproducible experiments."""
from .core import (Spectrum, CayleyGeometry, cayley_geometry, ordinal_probabilities,
                   spectra_from_probabilities, spectra_from_series, total_variation,
                   check_spectral_identities, keller_order_m3, path_energy)
from .graph import cayley_graph, plot_cayley_graph
__version__="1.0.0"
__all__=["Spectrum","CayleyGeometry","cayley_geometry","ordinal_probabilities",
         "spectra_from_probabilities","spectra_from_series","total_variation",
         "check_spectral_identities","keller_order_m3","path_energy",
         "cayley_graph","plot_cayley_graph"]
