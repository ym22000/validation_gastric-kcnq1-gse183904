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

## Readouts

The main readout is KCNQ1 detection rate:

- KCNQ1-positive cell: raw count `> 0`
- detection rate: percentage of KCNQ1-positive cells in each group

Expression intensity among KCNQ1-positive cells is reported as a secondary readout only. For sparsely detected genes such as ion channels, the fraction of positive cells is more informative than expression intensity among already positive cells.

## Statistical tests

- Detection rate: Fisher exact test on KCNQ1-positive versus KCNQ1-negative cells.
- Expression among KCNQ1-positive cells: two-sided Mann-Whitney/Wilcoxon test on `log1p(CPM)`.
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

The detection difference is supported by Fisher exact tests. Expression intensity among KCNQ1-positive cells is less informative and more similar between normal and intestinal tumor cells.

## Outputs

- `outputs/final/kcnq1_lauren_downsampled/KCNQ1_lauren_downsampled_analysis.pdf`
- `outputs/final/kcnq1_lauren_downsampled/KCNQ1_lauren_downsampled_analysis.xlsx`
- `outputs/final/kcnq1_lauren_downsampled/KCNQ1_lauren_downsampling_iteration_summary.tsv`
- `outputs/final/kcnq1_lauren_downsampled/KCNQ1_lauren_downsampled_representative_tests.tsv`

## Script

- `scripts/09_kcnq1_lauren_downsampling.py`
