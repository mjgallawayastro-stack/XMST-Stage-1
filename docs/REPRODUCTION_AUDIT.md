# XMST Paper I: 100-field reproduction audit

## Reproduction status

The controlled 100-field empirical-template Monte Carlo campaign completed successfully using:

- `quintana_2025_ob_stars_full(2).csv` (24,706 rows)
- `selected_template_groups.csv`
- a reconstructed `q25_group_membership.csv` generated directly from the nominal Q25 Delaunay–MST–Percolation–Jenks analysis.

The reconstructed membership is exact for the eight selected template labels:
192, 183, 109, 57, 51, 43, 37 and 25 stars (697 total).
The nominal fracture scale was 17.0483885 pc, consistent with the previously recorded 17.048476 pc difference being only rounding/numerical precision.

## Reproduced headline results

- Global detection of all eight templates: 100/100 fields.
- Global mean completeness: 0.9994.
- Global mean purity: 0.7445.
- Deliberate overlap resolved spatially: 2/100 fields.
- Deliberate overlap resolved after A_V pruning: 45/100 fields.
- A_V pruning resolved 43/98 initially merged overlaps.
- Final detection of all eight templates: 79/100 fields.
- Final mean completeness: 0.9454.
- Final mean purity: 0.8071.
- Median Jenks threshold: 17.028 pc.
- Mean Jenks threshold: 17.005 pc.
- Between-field Jenks standard deviation: 0.071 pc.
- Median within-field bootstrap standard deviation: 0.075 pc.
- Median percolation operational cut: 29.691 pc.
- Accepted A_V branch cuts: 177 in 78 fields.

These values match the previous archived report.

## Code-to-method audit

1. Each field contains 24,009 background stars plus 697 injected template stars, keeping the total at 24,706 rather than adding 697 stars to the full catalogue.
2. The eight templates are Q25 XMST structures extracted from the same parent catalogue, not independently labelled Q26 associations. The paper should not call this an independent Q26 ground-truth test.
3. The first two (largest) templates are placed on the same sightline with centre separation drawn uniformly from 57–83 pc.
4. The Delaunay graph, exact Euclidean MST, percolation limit and two-class Jenks threshold are recomputed independently for every outer field.
5. The Jenks uncertainty test uses 200 resamples with replacement of the subcritical edge-length distribution.
6. Retained components contain at least 10 stars.
7. Completeness is true recovered members / injected members.
8. Purity is true recovered members / total membership of the best-matching recovered group.
9. The controlled 100-field driver does not perturb distance uncertainties. Distance-uncertainty propagation is a separate nested analysis and has not been rerun in this reproduction package.
10. The revised non-Gaussian distance sampler should therefore be described and tested separately from these 100 baseline fields.

## Reproducibility note

The original membership table was not supplied, but it was reconstructed deterministically from the nominal catalogue using the same production graph and fracture-scale code. All eight selected component labels and sizes matched the template summary exactly.
