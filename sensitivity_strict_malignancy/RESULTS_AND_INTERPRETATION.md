# Results and interpretation

## Why this control was added

Diffuse and EMT programs can overlap with fibroblast and mesothelial transcriptional programs. The main Dataset 4 object was therefore kept unchanged and analyzed again under two stricter definitions. This tests whether the KCN result depends on the ambiguous clusters.

## Cohorts

### Original final

The original 6,573 high-confidence state cells were retained:

- 4,170 intestinal-like cells.
- 2,403 diffuse/EMT-like cells.
- 20 patients represented in both states.

### Cluster-screened

Clusters 2, 5 and 16 were removed because they combined low CNV support with mesothelial markers. Clusters 17 and 18 were removed because they combined low CNV support with strong fibroblast programs. These labels describe a conservative quality-control decision and do not prove that every excluded cell is non-malignant.

The resulting cohort contains:

- 4,167 intestinal-like cells.
- 631 diffuse/EMT-like cells.
- 20 patients represented.

### Orthogonal strict

A cell was retained when it came from a primary tumor, detected at least two epithelial-core genes, detected fewer than two mesothelial-core genes, and met one of the following conditions:

- `cnv_high = TRUE`; or
- tumor-program support together with an epithelial score higher than the fibroblast score, at least four epithelial-core genes, and fewer than three fibroblast-core genes.

The gene panels were:

- Epithelial core: `KRT8`, `KRT18`, `KRT19`, `EPCAM`, `CDH1`, `MUC1`, `TACSTD2`.
- Fibroblast core: `COL1A1`, `COL1A2`, `COL3A1`, `DCN`, `LUM`, `FN1`, `SPARC`.
- Mesothelial core: `CALB2`, `WT1`, `UPK3B`, `LRRN4`, `MSLN`.

This conservative cohort contains:

- 4,000 intestinal-like cells.
- 419 diffuse/EMT-like cells.
- 20 patients represented.

## KCNQ1

`KCNQ1` remains preferentially detected in intestinal-like cells under every definition:

| Cohort | Intestinal detection | Diffuse/EMT detection | Fisher OR | FDR |
| --- | ---: | ---: | ---: | ---: |
| Original final | 25.7% | 1.7% | 20.39 | 7.07e-178 |
| Cluster-screened | 25.7% | 5.2% | 6.25 | 3.12e-37 |
| Orthogonal strict | 25.6% | 4.3% | 7.64 | 1.45e-28 |

The patient-level comparison has the same direction. In the orthogonal-strict cohort, 8 of 10 eligible patients show higher intestinal detection, one shows the opposite direction, and one is tied. The paired Wilcoxon FDR is 0.0193.

The within-patient rank association also remains positive with the intestinal score and negative with the diffuse/EMT score. The correlations are weak, but their direction is stable across all cohorts.

## KCNE3

`KCNE3` also remains preferentially detected in intestinal-like cells:

| Cohort | Intestinal detection | Diffuse/EMT detection | Fisher OR | FDR |
| --- | ---: | ---: | ---: | ---: |
| Original final | 28.2% | 5.7% | 6.43 | 1.35e-122 |
| Cluster-screened | 28.1% | 10.1% | 3.47 | 3.33e-25 |
| Orthogonal strict | 27.8% | 11.0% | 3.11 | 1.99e-15 |

All 10 eligible patients in the orthogonal-strict cohort show higher intestinal detection. The paired Wilcoxon FDR is 0.0178. Continuous-score associations remain positive for the intestinal program and negative for the diffuse/EMT program within patients.

## KCNE2

`KCNE2` is detected more frequently in intestinal-like cells in all three cohorts, but the signal is small:

- Detection differences range from 1.4 to 3.3 percentage points.
- The cluster-screened Fisher FDR is 0.099.
- The within-patient association with the intestinal score is weak and slightly negative.

The evidence is therefore insufficient to describe `KCNE2` as robustly associated with the intestinal state in this dataset.

## Robustness assessment

The main biological direction is robust for `KCNQ1` and `KCNE3` because it is supported by four complementary observations:

1. Detection is higher in intestinal-like cells in the original and both stricter cohorts.
2. Fisher odds ratios remain greater than 1 after ambiguous clusters are removed.
3. Most eligible patients show the same direction independently.
4. Within-patient correlations remain positive for the intestinal score and negative for the diffuse/EMT score.

The reduction in odds ratios after screening shows that ambiguous low-KCN cells increased the original contrast. Therefore, the original effect size should not be interpreted literally. The defensible conclusion is that `KCNQ1` and `KCNE3` are preferentially retained or detected in intestinal-like malignant epithelial states, not that they are exclusive markers of those states.

The strict analysis favors specificity and may remove genuine EMT tumor cells with weak CNV or reduced epithelial markers. It is therefore a sensitivity control rather than a replacement for the biologically heterogeneous main object.

