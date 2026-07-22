# KCNQ1 Lauren-level downsampling analysis

## Question

Does KCNQ1 detection differ between epithelial cells from primary normal gastric samples, intestinal-type gastric tumors and diffuse-type gastric tumors?

## Data used

The analysis uses epithelial candidate cells from GSE183904.

The three groups are:

- `Primary normal epithelial`: epithelial candidate cells from `Primary_Normal` samples.
- `Tumor intestinal patient`: epithelial candidate cells from `Primary_Tumor` samples annotated as intestinal by Lauren type.
- `Tumor diffuse patient`: epithelial candidate cells from `Primary_Tumor` samples annotated as diffuse by Lauren type.

This is a Lauren sample-level analysis. It does not classify each tumor cell as intestinal-like or diffuse/EMT-like by transcriptional signatures.

## Why downsampling

The three groups contain different numbers of cells. To make the comparison more balanced, each group was randomly downsampled to the size of the smallest group.

- smallest group: `Tumor diffuse patient`
- downsample size: `3,368` cells per group
- iterations: `500`
- random seed: `12345`

Downsampling balances the number of displayed cells, but does not make cells from the same patient independent. The pooled-cell tests are therefore exploratory only.

## Readouts

The main readout is KCNQ1 detection rate:

- KCNQ1-positive cell: raw count `> 0`
- detection rate: percentage of KCNQ1-positive cells in each group

Expression intensity among KCNQ1-positive cells is reported as a secondary readout only. For sparsely detected genes such as ion channels, the fraction of positive cells is more informative than expression intensity among already positive cells.

## Statistical tests

- Detection rate: Fisher exact test on KCNQ1-positive versus KCNQ1-negative cells.
- Normal versus tumor inference: exact paired Wilcoxon signed-rank test across the five patients with matched primary normal and primary tumor samples.
- Intestinal versus diffuse inference: two-sided Mann-Whitney test on patient-level summaries (`14` intestinal and `6` diffuse patients).
- Expression among KCNQ1-positive cells: exploratory two-sided Mann-Whitney test on `log1p(CP10K)`, where CP10K is the raw count divided by the cell library size and multiplied by `10,000`.
- Multiple testing: Benjamini-Hochberg FDR across pairwise comparisons.

## Main result

After downsampling, KCNQ1 detection follows the expected direction:

```text
Primary normal epithelial > Tumor intestinal patient > Tumor diffuse patient
```

Representative downsample:

| Group | Cells | KCNQ1-positive cells | Detection |
|---|---:|---:|---:|
| Primary normal epithelial | 3,368 | 480 | 14.25% |
| Tumor intestinal patient | 3,368 | 370 | 10.99% |
| Tumor diffuse patient | 3,368 | 262 | 7.78% |

Pooled-cell Fisher tests support the displayed direction, but they are not the main inference because cells are nested within patients. The matched-patient and patient-level tests are the appropriate inferential results. Expression intensity among KCNQ1-positive cells is secondary and more similar between normal and intestinal tumor cells.

## Matched-patient result

Five patients have both a primary normal epithelial sample and a primary tumor sample. These pairs are tested with an exact paired Wilcoxon signed-rank test. The small number of pairs limits statistical power, so effect direction and paired differences are reported together with the p-value.

KCNQ1 detection was higher in the normal sample for four of the five pairs. The median paired difference was `+6.89` percentage points, but the exact paired test was not significant (`p = 0.625`). Mean `log1p(CP10K)` showed the same direction in four pairs and was also not significant (`p = 0.3125`).

## Outputs

- `outputs/final/kcnq1_lauren_downsampled/KCNQ1_lauren_downsampled_analysis.pdf`
- `outputs/final/kcnq1_lauren_downsampled/KCNQ1_lauren_downsampled_analysis.xlsx`
- `outputs/final/kcnq1_lauren_downsampled/KCNQ1_lauren_downsampling_iteration_summary.tsv`
- `outputs/final/kcnq1_lauren_downsampled/KCNQ1_lauren_downsampled_representative_tests.tsv`
- `outputs/final/kcnq1_primary_normal_vs_tumor/kcnq1_matched_normal_tumor_pairs.tsv`
- `outputs/final/kcnq1_primary_normal_vs_tumor/kcnq1_matched_normal_tumor_paired_tests.tsv`
- `outputs/final/kcnq1_primary_normal_vs_tumor/kcnq1_tumor_lauren_patient_tests.tsv`

## Script

- `scripts/09_kcnq1_lauren_downsampling.py`
