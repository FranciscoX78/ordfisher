# Publishing version 1.0.0: checklist

**Do not announce a DOI, a PyPI release, or a public GitHub repository before these events actually occur.**

1. **Scientific sign-off:** confirm that the paper's current approved figures and numerical descriptions align with the included scripts and archived revision. Preserve the exact certificate provenance distinction.
2. **Credit and governance:** the named authors are Giovanni Campanini, Fernando Pose, Carlos Garcia, and Francisco O. Redelico; Francisco O. Redelico is the maintainer. Confirm the approval of the license and copyright notice by all relevant rightsholders before release.
3. **Final QA:** run `python -m pytest -q`, the exact Markov verifier, and AR(1) verifier. Run the `ordpy` comparison on an environment where `ordpy` is installed. Test supported Python versions via continuous integration. Run the complete logistic/Markov pipelines if the deposit is intended to include their fresh CSV outputs.
4. **Build:** use `python -m pip install build twine`, `python -m build`, and `python -m twine check dist/*` in an environment with those packaging tools available. Test installation from the wheel in a fresh environment.
5. **GitHub:** create a public repository, commit the reviewed source (no `results/` folders or untracked private data), run CI, and tag the reviewed commit `v1.0.0`. Confirm the repository's final address before adding it to metadata.
6. **Zenodo:** connect Zenodo to the actual repository, archive the **specific v1.0.0 release** and obtain its version-specific DOI. Enter the real DOI and URL in `CITATION.cff` / release documentation, then update the paper's code availability statement. The archive should identify the code version used for the manuscript.
7. **PyPI (optional):** reserve/confirm the distribution name and publish only after checking package availability. A Python import name does not guarantee that the distribution name on PyPI is free. Do not claim `pip install ordfisher` works until tested after publication.
8. **Manuscript:** cite the **archived release DOI**, not a fabricated placeholder, and clarify the new seeds and newly reconstructed Markov certificate.

## Suggested concise code-availability statement (after a real DOI exists)

> The `ordfisher` Python library and the scripts used to reproduce the spectral calculations and numerical examples are openly available from the project's public repository, with the version used in this work archived under the cited version-specific DOI. The supplied rational Markov certificate is a newly constructed exact realization satisfying the stated constraints; simulations use recorded random seeds.

Insert the *actual* links and DOI only after publication.
