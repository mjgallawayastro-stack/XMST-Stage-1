# XMST nested distance-uncertainty experiment

Completed outer fields: **100**
Uncertainty realisations per field: **100**
Total perturbed XMST reconstructions: **10000**

## Baseline reproduction
- Median nominal Jenks scale: 17.027820 pc
- Between-field nominal Jenks SD: 0.071433 pc
- Mean nominal template completeness: 0.999441
- Mean nominal template purity: 0.744453

## Distance-uncertainty result
- Median perturbed Jenks scale: 17.250891 pc
- Perturbed Jenks SD: 0.304939 pc
- Median Jenks shift relative to each nominal field: 0.245293 pc
- Median absolute Jenks shift: 0.282423 pc
- Mean nominal-group Jaccard: 0.433191
- Median group median-Jaccard: 0.444664
- Mean group recovery fraction (J >= 0.50): 0.500243
- Mean group recovery fraction (J >= 0.75): 0.246505
- Mean group recovery fraction (J >= 0.90): 0.072216
- Median star membership persistence: 0.740000
- Stars with persistence >= 0.80: 0.447479
- Mean injected-template completeness under uncertainty: 0.860530
- Mean injected-template purity under uncertainty: 0.712601
- Fraction of template cases satisfying >=50% detection: 0.992975

## Sampling definition
For every synthetic star, the lower and upper *fractional* Q25 distance uncertainties are preserved.
The moved star's nominal MC radius is treated as its median distance, and synthetic 16th/84th percentiles
are constructed by applying the source star's fractional lower/upper uncertainty. Distances are then drawn
from the same percentile-matched two-sided normal approximation used in the Q25 stability code.
