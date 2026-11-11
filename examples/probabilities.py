"""Calculate OPS and OFS from a known lexicographically ordered S3 law."""
from ordfisher import spectra_from_probabilities

p = [1/3, 1/15, 2/15, 1/5, 4/15, 0]
r = spectra_from_probabilities(p, m=3)
print('frequencies:', r.frequencies)
print('OPS:', r.ops)
print('OFS:', r.ofs)
print('Fisher energy:', r.fisher_energy)
