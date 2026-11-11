# Reproducing the accompanying OFS paper

This repository provides code for both the reusable OPS/OFS library and the companion **research examples**. The examples do not alter the core spectral definitions.

## 1. Fast validation (seconds to minutes)

```bash
python -m pytest -q
python scripts/verify_markov.py --out results/markov/verification.json
python scripts/verify_ar1.py --out results/ar1
python scripts/run_logistic.py --out results/logistic_quick --seed 7 \\
  --r-min 3.70 --r-max 3.71 --r-step 0.01 \\
  --repetitions 2 --burn 200 --retained 500
```

The short logistic example is a **software smoke test**, not evidence of reproducing the entire paper's statistical results.

## 2. Full logistic-map protocol

```bash
python scripts/run_logistic.py --out results/logistic --seed 20261105
```

This samples `r=3.50,...,4.00` in steps of `0.01`: 51 parameter values, ten initial conditions each, 5,000 discarded transients and 120,000 retained observations per trajectory, with `m=4` and `tau=1`. A total of 510 trajectories are evaluated. The script writes individual CSV rows, aggregated means and sample standard deviations, metadata, and plots.

**Important:** the seed is newly assigned and recorded. Historical initial conditions were not recovered; the experiment is repeatable from this archive, but is not expected to reproduce the *exact original* simulated samples bit for bit.

## 3. Exact Markov certificate

```bash
python scripts/verify_markov.py --out results/markov/verification.json
```

The bundled `ordfisher/data/markov_q_exact.json` contains *two newly reconstructed* strictly positive rational four-block laws (`8^4 = 4096` numerators each, with a common denominator per process). Integer-arithmetic checks cover normalization, uniform incoming/outgoing three-block marginals, the required lag-three covariance identity, and all 24 induced ordinal probabilities. Algebraic graph arguments verify matching OPS and scalar Fisher energies; a numerical spectral calculation verifies nonzero OFS separation.

It is a **new certificate for the manuscript's exact specified ordinal distributions**, **not** the original unpublished certificate file. The optional `scripts/rebuild_rational_certificate.py` can generate another valid certificate after installing SymPy; the distributed one is the reproducible fixed artifact.

## 4. Stationary Markov simulations

```bash
python scripts/run_markov_simulation.py --out results/markov --seed 20261105
```

The simulation uses 30 realizations per process and 32,768 overlapping ordinal windows per realization. Sample statistics depend on the chosen seed and newly reconstructed rational chains; they should be compared with the exact stationary probabilities rather than forced to match historical realization-level values.

## 5. Exact Gaussian AR(1) illustration

```bash
python scripts/verify_ar1.py --out results/ar1
```

For `m=4`, AR(1) pattern probabilities are calculated using a Gaussian orthant identity derived from the covariance function. The script computes the `phi=-0.5` and positive-`phi` pair with equal graph Fisher energy and a parameter sweep; neither uses Monte Carlo estimates.

## 6. Figures and source assets

- `plot_cayley_graph(4, layout='appendix')` regenerates the S4 Cayley figure from frozen original placement and styling data.
- The experimental scripts generate only their named example figures and the corresponding numerical tables; no all-purpose diagnostic figure generator is included in the v1.0.0 public package.
- The approved Figure 1 of the paper is a historical illustration, not a required output of the main library. Generated redraws of that illustration are not asserted to be identical artwork.

## Reference outcomes (not a substitute for running the code)

| Quantity | Expected approximate result |
|---|---:|
| Markov normalized OFS total-variation separation | `0.23122927995` |
| Common Markov Cayley Fisher energy | `0.06155206485` |
| Positive AR(1) parameter with same Fisher energy as `-0.5` | `0.6536965296` |
| AR(1) Fisher-energy allocation distance | `0.4127584220` |
| Logistic OPS/OFS distance around `r=3.70` for the recorded new seed | `0.13948` |
| Logistic OPS/OFS distance around `r=3.93` for the recorded new seed | `0.20500` |

These numbers do not prove general classification advantage or replace a careful study of finite-sample uncertainty; the ordinal windows overlap.
