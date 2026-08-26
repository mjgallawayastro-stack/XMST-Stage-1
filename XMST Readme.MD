{\rtf1\ansi\ansicpg1252\cocoartf2870
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 # XMST\
\
**XMST: An Extended Minimum Spanning Tree Framework with Objective Fracture-Scale Selection**\
\
XMST is a spatial clustering framework designed to identify candidate stellar structures using a global minimum spanning tree (MST) together with an objective, data-driven fracture-scale criterion.\
\
The method combines **percolation analysis** with **Jenks Natural Breaks optimisation**. Percolation is first used to restrict the MST edge-length distribution to the subcritical regime, after which Jenks optimisation determines the final fracture scale.\
\
This repository contains the code, configuration files, validation outputs and manuscript material associated with the XMST methodology paper.\
\
**Author:** Mark J. Gallaway\
\
---\
\
## Paper\
\
This repository accompanies:\
\
**Gallaway, M. J. (2026),  \
\'93XMST: An Extended Minimum Spanning Tree Framework with Objective Fracture-Scale Selection\'94**\
\
Submitted manuscript details and DOI will be added following publication.\
\
The version of the software used for the paper will also be permanently archived on Zenodo.\
\
---\
\
## Method overview\
\
The principal XMST workflow is:\
\
1. Construct a three-dimensional representation of the stellar catalogue.\
2. Generate a three-dimensional Delaunay triangulation.\
3. Recover the global Euclidean minimum spanning tree.\
4. Perform percolation analysis on the ordered MST edge lengths.\
5. Identify the transition towards global connectivity.\
6. Restrict the edge-length distribution to the subcritical regime.\
7. Apply two-class Jenks Natural Breaks optimisation.\
8. Adopt the resulting boundary as the XMST fracture scale.\
9. Fracture the global MST.\
10. Retain subtrees satisfying the adopted minimum membership criterion.\
11. Assess recovery and membership stability through Monte Carlo and distance-uncertainty experiments.\
\
Additional reddening-based diagnostics were investigated as a post-processing method for identifying possible line-of-sight superpositions.\
\
---\
\
## Validation\
\
XMST was tested using **1000 controlled Monte Carlo realisations**, each containing eight empirical stellar-association templates embedded within a realistic stellar background.\
\
The principal validation results were:\
\
- median Percolation\'96Jenks fracture scale: **17.028 pc**\
- between-realisation standard deviation: **0.066 pc**\
- mean injected-structure completeness: **0.9996**\
- mean purity: **0.7413**\
- all eight injected structures satisfied the main detection criterion in all 1000 nominal realisations\
\
Two templates were deliberately placed in a difficult overlapping configuration. Spatial XMST alone resolved this pair in **15/1000** realisations. A reddening-based post-processing diagnostic increased the number of resolved cases to **446/1000**, although with a corresponding completeness\'96purity trade-off.\
\
The published distance uncertainties were also propagated through **10,000 additional XMST reconstructions**. The median absolute change in fracture scale was **0.282 pc**, while **99.30%** of injected-template cases continued to satisfy the adopted detection criterion.\
\
The repository also contains direct comparisons with the traditional **CDF** and **Mean Edge** fracture-scale criteria.\
\
---\
\
## Repository structure\
\
```text\
XMST/\
\uc0\u9500 \u9472 \u9472  code/\
\uc0\u9474    \u9500 \u9472 \u9472  percolation_jenks/\
\uc0\u9474    \u9500 \u9472 \u9472  cdf/\
\uc0\u9474    \u9500 \u9472 \u9472  mean_edge/\
\uc0\u9474    \u9492 \u9472 \u9472  distance_uncertainty/\
\uc0\u9474 \
\uc0\u9500 \u9472 \u9472  results/\
\uc0\u9474    \u9500 \u9472 \u9472  percolation_jenks/\
\uc0\u9474    \u9500 \u9472 \u9472  cdf/\
\uc0\u9474    \u9500 \u9472 \u9472  mean_edge/\
\uc0\u9474    \u9492 \u9472 \u9472  distance_uncertainty/\
\uc0\u9474 \
\uc0\u9500 \u9472 \u9472  paper/\
\uc0\u9474    \u9500 \u9472 \u9472  Paper1.tex\
\uc0\u9474    \u9500 \u9472 \u9472  references.bib\
\uc0\u9474    \u9492 \u9472 \u9472  figures/\
\uc0\u9474 \
\uc0\u9500 \u9472 \u9472  data/\
\uc0\u9474    \u9492 \u9472 \u9472  README.md\
\uc0\u9474 \
\uc0\u9500 \u9472 \u9472  docs/\
\uc0\u9474    \u9500 \u9472 \u9472  REPRODUCTION_AUDIT.md\
\uc0\u9474    \u9500 \u9472 \u9472  RELEASE_CHECKLIST.md\
\uc0\u9474    \u9492 \u9472 \u9472  ZENODO_DATA_MANIFEST.md\
\uc0\u9474 \
\uc0\u9500 \u9472 \u9472  requirements.txt\
\uc0\u9500 \u9472 \u9472  SHA256SUMS.txt\
\uc0\u9492 \u9472 \u9472  README.md}