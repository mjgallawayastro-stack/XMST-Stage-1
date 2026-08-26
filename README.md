# XMST — Paper I reproducibility package

This repository candidate contains the code, configuration files, compact validation outputs, and manuscript material associated with **XMST: An Extended Minimum Spanning Tree Framework with Objective Fracture-Scale Selection**.

## Status

This is a pre-publication staging package. It should be checked before the repository is made public or archived.

## Analyses included

- Percolation–Jenks XMST validation over 1,000 controlled Monte Carlo fields.
- CDF fracture-scale comparison over the same 1,000 fields.
- Mean Edge fracture-scale comparison over the same 1,000 fields.
- Nested distance-uncertainty analysis: 100 nominal fields × 100 perturbations = 10,000 XMST reconstructions.

## Software requirements

Python 3.10+ is recommended. The production distance-uncertainty run recorded Python 3.14.5, NumPy 2.4.6, pandas 3.0.5 and SciPy 1.18.0 on macOS/Apple Silicon. The Monte Carlo comparison code also uses scikit-learn and matplotlib.

Install the core dependencies with:

```bash
python3 -m pip install -r requirements.txt
```

## Repository layout

- `code/` — production and comparison scripts.
- `results/` — compact headline outputs and method configuration files.
- `paper/` — manuscript source, bibliography and final paper figures.
- `docs/` — reproduction/audit notes.
- `data/` — instructions for obtaining or staging the Q25 input catalogue.

## Important data note

The full Q25-derived input catalogue is deliberately **not** included in this GitHub candidate until its redistribution terms have been checked. The original working archive contains the exact input files used by the production runs. Full Monte Carlo output tables are intended for the permanent Zenodo archive rather than the GitHub repository.

## Reproduction

The original production runners expect their input CSV files in an `input/` directory alongside the relevant scripts. See `data/README.md` and the analysis-specific code for the expected filenames.

## Citation and DOI

A `CITATION.cff` file and Zenodo DOI should be added once the repository URL and archived release exist.
