# Chronological analysis pipeline

## 1. Data scope

The input is the author-processed raw gene-count archive from GSE183904. The published Cell Ranger and quality-control criteria were 500 to 6,000 detected genes per cell and at most 20% mitochondrial reads. The GEO matrices already passed this author-level filtering.

The present analysis retains 20 primary tumors with an available Lauren diagnosis: 14 intestinal and 6 diffuse tumors. Five matched primary normal samples are added only to construct the inferCNV reference. Mixed tumors, metastases, peritoneal samples and samples without a usable Lauren label are excluded.

| Step | Cells |
|---|---:|
| Cells in the selected GEO matrices | 113,470 |
| Broad epithelial candidates | 28,845 |
| Malignant epithelial candidates | 11,328 |
| After immune-contamination filtering | 8,901 |
| Final high-confidence states | 6,573 |

## 2. Broad epithelial candidate gate

**What is done.** The compressed count matrices are streamed sample by sample. A cell is retained when at least three epithelial core genes are detected, or when at least two are detected and the epithelial score is greater than 65% of the strongest immune, fibroblast, endothelial or pericyte score.

**Why.** A strict EPCAM-only gate would lose plastic or partially mesenchymal tumor cells. The relaxed branch keeps possible EMT cells while still controlling obvious non-epithelial profiles.

**Core genes.** EPCAM, KRT7, KRT8, KRT18, KRT19, KRT20, MUC1, CDH1 and TACSTD2.

**Output.** 28,845 broad epithelial candidates from 113,470 cells. The sparse count matrix is written in Matrix Market format to avoid loading the complete atlas as a dense matrix.

## 3. Seurat representation

**What is done.** Counts are normalized with `LogNormalize` and a scale factor of 10,000. Two thousand variable genes are selected. Data are scaled while regressing total UMI counts. Thirty principal components are computed; PCs 1 to 20 are used for the neighbor graph and t-SNE. The graph uses 30 neighbors and clustering resolution 0.8. The random seed is 12345.

**Why.** This representation summarizes transcriptional variation among epithelial candidates and provides a reproducible visualization. The t-SNE is descriptive and is not used alone to define malignancy or intestinal/diffuse identity.

## 4. RNA-inferred CNA support

**What is done.** inferCNV compares primary-tumor epithelial candidates with 750 matched-normal epithelial reference cells, balanced at 150 cells from each of five normal patients. Genes detected in fewer than 20 cells are removed before inferCNV. The run uses the GRCh38 gene order, cutoff 0.1, denoising, no HMM and eight threads.

**Threshold.** A CNA score is considered high when it exceeds the 95th percentile of the normal reference distribution: 0.0552533.

**Why.** Tumor cells can lose epithelial markers during plasticity. RNA-inferred CNA provides an orthogonal rescue signal, but is not used by itself because transcriptional stress and technical effects can also affect inferred CNA profiles.

**Difference from the source article.** Kumar et al. used CONICSmat together with epithelial lineage information. This project uses inferCNV, as in the dataset2 workflow, because the exact author cell-level annotation and analysis code are not provided in GEO. It is therefore an adapted reconstruction, not an exact reproduction.

## 5. Malignant epithelial selection

**What is done.** Malignant and non-malignant epithelial programs are scored from Zhou et al. Supplementary Table S4. Genes are retained from that table when adjusted p-value is at most 0.01, average log fold-change is at least 0.5 and detection is at least 25%.

The primary selection requires a tumor sample and a tumor-program difference above the 95th percentile measured in normal epithelial cells (1.10733). Additional cells can be rescued when their inferCNV score is above the normal 95th percentile. Cells with at least two coherent immune markers are removed.

**Selection result.** 8,543 cells are selected by the tumor program and 2,785 additional cells are rescued by inferCNV. After removing 2,427 immune-contaminated observations, 8,901 malignant epithelial candidates remain.

## 6. Intestinal and diffuse/EMT state scores

**What is done.** UCell scores are calculated with a maximum rank of 1,500. The intestinal program contains 22 genes. The diffuse/EMT score combines a 36-gene EMT-up program with loss of a 10-gene epithelial-junction program:

```text
Diffuse/EMT score = z(EMT-up UCell score) - z(epithelial-junction UCell score)
State delta       = z(Intestinal score) - z(Diffuse/EMT score)
```

KCNQ1, KCNE2 and KCNE3 are absent from all classifiers, preventing circular classification.

**High-confidence rule.** A cell is labeled only when the absolute state delta is at least 0.5 standard deviations and the dominant score is not below zero. Other cells remain indeterminate.

**Result.** 4,170 intestinal-like cells, 2,403 diffuse/EMT-like cells and 2,328 indeterminate cells. The retained high-confidence fraction is 73.85% of malignant epithelial candidates.

## 7. KCN detection analysis

**What is measured.** Detection is defined as at least one raw transcript count in a cell. The workbook reports expressing cells, total cells, raw transcript sums and mean raw counts.

**Tests.** A two-sided Fisher exact test compares expressing and non-expressing cells between intestinal-like and diffuse/EMT-like states. Benjamini-Hochberg correction is applied across the three KCN genes.

**Why.** Fisher exact test is a conventional bioinformatics choice for comparing detection rates in sparse single-cell data, because the readout is binary at cell level: detected or not detected. Patient-level heterogeneity is then checked with the continuous within-patient rank analysis below.

## 8. Signature correlations

**What is measured.** Spearman correlation relates normalized KCN expression to the continuous intestinal and diffuse/EMT scores.

**Tests.** Pooled Spearman correlations are descriptive. The patient-aware result centers ranks within each patient and estimates significance with 1,000 stratified permutations using seed 12345. Benjamini-Hochberg correction is applied across gene-signature pairs.

**Why.** Spearman is suitable for sparse, non-Gaussian single-cell values and tests monotonic association without assuming linearity. The within-patient analysis asks whether the direction persists beyond patient composition.

## 9. Independent clinical sanity check

Lauren diagnosis is not used to classify individual cells. After classification, the diffuse/EMT-like fraction is compared across patient diagnoses. The median is 0.670 in diffuse patients and 0.289 in intestinal patients; Wilcoxon p = 0.0433. This supports biological coherence, while the patient-level heterogeneity remains visible.
