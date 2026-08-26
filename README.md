# XMST

**XMST: An Extended Minimum Spanning Tree Framework with Objective Fracture-Scale Selection**

XMST is a spatial clustering framework for identifying candidate stellar structures from a global minimum spanning tree (MST) using an objective, data-driven fracture-scale criterion.

The method combines **percolation analysis** with **Jenks Natural Breaks optimisation**. Percolation first restricts the MST edge-length distribution to the subcritical regime; two-class Jenks optimisation then determines the final XMST fracture scale.

This repository contains the code, configuration files, compact validation outputs, and manuscript material associated with the XMST methodology paper.

**Author:** Mark J. Gallaway

## Paper

This repository accompanies:

**Gallaway, M. J. (2026), “XMST: An Extended Minimum Spanning Tree Framework with Objective Fracture-Scale Selection”**

Publication details and DOI will be added when available. A frozen release of the code and full validation outputs will also be archived on Zenodo.

## Method overview

The principal XMST workflow is:

1. Represent the stellar catalogue in three-dimensional Cartesian space.
2. Construct a three-dimensional Delaunay triangulation.
3. Recover the global Euclidean minimum spanning tree.
4. Perform percolation analysis on the ordered MST edge lengths.
5. Identify the transition towards global connectivity.
6. Restrict the edge-length distribution to the subcritical regime.
7. Apply two-class Jenks Natural Breaks optimisation.
8. Adopt the resulting boundary as the XMST fracture scale.
9. Fracture the global MST.
10. Retain subtrees satisfying the adopted minimum membership criterion.
11. Assess recovery and membership stability through Monte Carlo and distance-uncertainty experiments.

Reddening-based diagnostics were additionally investigated as a post-processing method for identifying possible line-of-sight superpositions.

## Validation

XMST was tested using **1000 controlled Monte Carlo realisations**, each containing eight empirical stellar-association templates embedded within a realistic stellar background.

Headline results:

- median Percolation–Jenks fracture scale: **17.028 pc**
- between-realisation standard deviation: **0.066 pc**
- mean injected-structure completeness: **0.9996**
- mean purity: **0.7413**
- all eight injected structures satisfied the main detection criterion in all 1000 nominal realisations

Two templates were deliberately placed in a difficult overlapping configuration. Spatial XMST alone resolved this pair in **15/1000** realisations. A reddening-based post-processing diagnostic increased the number of resolved cases to **446/1000**, with a corresponding completeness–purity trade-off.

Published stellar-distance uncertainties were propagated through **10,000 additional XMST reconstructions**. The median absolute change in fracture scale was **0.282 pc**, while **99.30 per cent** of injected-template cases continued to satisfy the adopted detection criterion.

The repository also contains direct comparisons with the traditional **CDF** and **Mean Edge** fracture-scale criteria.

## Repository structure

```text
XMST-Stage-1/
├── code/
│   ├── percolation_jenks/
│   ├── cdf/
│   ├── mean_edge/
│   └── distance_uncertainty/
├── results/
│   ├── percolation_jenks/
│   ├── cdf/
│   ├── mean_edge/
│   └── distance_uncertainty/
├── paper/
│   ├── Paper1.tex
│   ├── references.bib
│   └── figures/
├── data/
│   └── README.md
├── docs/
│   ├── REPRODUCTION_AUDIT.md
│   ├── RELEASE_CHECKLIST.md
│   └── ZENODO_DATA_MANIFEST.md
├── requirements.txt
├── SHA256SUMS.txt
├── LICENSE
└── README.md
```

### Code

`code/percolation_jenks/` contains the production Percolation–Jenks Monte Carlo validation code. Equivalent CDF and Mean Edge comparison experiments are contained in their respective directories. The nested distance-uncertainty analysis is in `code/distance_uncertainty/`.

### Results

`results/` contains compact summary outputs and configuration files used to verify the numerical results reported in the paper. Full Monte Carlo output tables are intended for the associated Zenodo archive rather than the GitHub repository.

## Installation

XMST was developed and tested using Python 3. The core dependencies are listed in `requirements.txt`.

```bash
python3 -m pip install -r requirements.txt
```

## Input data

The parent stellar sample used for Paper I is the catalogue of candidate O- and B-type stars presented by **Quintana, Wright & Martínez García (2025)**.

The production analyses expect the required input CSV files in an `input/` directory. See `data/README.md` for filenames and provenance.

The Q25 catalogue and derived input tables are not redistributed in this repository pending confirmation of their redistribution terms. Where necessary, users should obtain the source catalogue from the original publication and follow the accompanying documentation to reconstruct the inputs.

## Reproducibility

The Monte Carlo validation uses controlled random seeds and fixed analysis parameters. Compact production configurations, summary outputs, and SHA-256 checksums are provided in this repository. Additional information on the verified production runs is provided in `docs/REPRODUCTION_AUDIT.md`.

## Interpretation

XMST identifies **candidate stellar structures**, not definitive astrophysical memberships. The validation shows that the global fracture scale and the existence of recovered structures are highly stable, while individual membership — particularly in peripheral, locally sparse, or strongly overlapping regions — can be more sensitive to observational uncertainty.

Additional astrophysical information such as proper motions, radial velocities, photometry, spectroscopy, and stellar ages may therefore be required to refine individual memberships.

## Citation

If you use XMST in published work, please cite the associated paper. A `CITATION.cff` file and permanent Zenodo DOI will be added to the frozen Paper I software release.

## Licence

XMST is released under the MIT License. See `LICENSE`.

## Contact

Questions, bug reports, and reproducibility queries can be raised through the GitHub Issues page for this repository.
