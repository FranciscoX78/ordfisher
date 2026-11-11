# External concordance of ordinal distributions (`ordpy`)

**Date:** 2026-11-11  
**Environment reported by Google Colab:** Python 3.13.18, `ordfisher` 1.0.0, `ordpy` 1.2.3.  
**Status:** PASS (user-provided Colab console output; no machine-readable JSON from the Colab run was separately archived at preparation time).

The two libraries were given matching time series and parameters, and the resulting ordinal-pattern distributions were compared for the same lexicographically ordered patterns.

| Input | Pattern lengths | Delays | Results |
|---|---|---|---|
| Gaussian noise | 3, 4, 5 | 1, 2 | 6 comparisons; maximum absolute error 0 |
| Gaussian AR(1) | 3, 4, 5 | 1, 2 | 6 comparisons; maximum absolute error 0 |
| Logistic map | 3, 4, 5 | 1, 2 | 6 comparisons; maximum absolute error 0 |

The Colab run also printed **`PASS: 18 ordinal-distribution comparisons and m=4 spectral consistency.`**

This is useful independent validation of the ordinal-probability implementation under these tested conditions. It does **not** by itself independently validate the entire graph-spectral theory or all simulations. Run `python scripts/check_ordpy_agreement.py` with the optional `ordpy` dependency installed for an additional command-line smoke test.
