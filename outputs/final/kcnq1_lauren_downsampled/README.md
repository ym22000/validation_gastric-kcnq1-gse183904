# KCNQ1 Lauren downsampling output

This folder contains the final Lauren-level KCNQ1 comparison.

## Groups

- `Primary normal epithelial`: epithelial candidate cells from primary normal gastric samples.
- `Tumor intestinal patient`: epithelial candidate cells from intestinal-type primary tumors.
- `Tumor diffuse patient`: epithelial candidate cells from diffuse-type primary tumors.

## Method

- Downsampling to `3,368` cells per group.
- `500` random iterations.
- Seed: `12345`.
- Main metric: KCNQ1 detection rate, defined as raw count `> 0`.
- The pooled-cell Fisher test is exploratory because cells from one patient are not independent replicates.
- The main normal-versus-tumor inference is an exact paired Wilcoxon test across the five matched patients.
- Intestinal-versus-diffuse inference uses a patient-level Mann-Whitney test.
- Secondary metric: expression among KCNQ1-positive cells only, using `log1p(CP10K)`.

## Files

- `KCNQ1_lauren_downsampled_analysis.pdf`: final figure.
- `KCNQ1_lauren_downsampled_analysis.xlsx`: Excel workbook with method, summaries, tests and iterations.
- `KCNQ1_lauren_full_group_summary.tsv`: full non-downsampled group summary.
- `KCNQ1_lauren_downsampled_representative_summary.tsv`: one representative downsample summary.
- `KCNQ1_lauren_downsampled_representative_tests.tsv`: pairwise tests on the representative downsample.
- `KCNQ1_lauren_downsampling_iteration_summary.tsv`: median and 95% interval across 500 downsampling iterations.
- `KCNQ1_lauren_downsampling_iterations.tsv`: all iteration-level values.
