# Malignant-cell sensitivity check

The main Dataset 4 analysis contains heterogeneous intestinal-like and diffuse/EMT-like malignant epithelial states. EMT cells can share transcriptional programs with fibroblast or mesothelial cells. This additional control asks a simple question: do the KCN results remain the same when malignancy is defined more strictly?

The original object and the main pipeline were not changed.

## What was compared

Three versions of the final population were analyzed:

1. `Original final`: the 6,573 cells from the main analysis.
2. `Cluster-screened`: the same cells without clusters 2, 5, 16, 17 and 18. These clusters had weak CNV support and strong mesothelial or fibroblast-like programs.
3. `Orthogonal strict`: cells supported by CNV, or by a tumor program together with a dominant epithelial identity.

The stricter groups are sensitivity controls. They do not replace the original annotation because genuine EMT tumor cells can lose epithelial markers or show weak inferred CNV.

## Result

`KCNQ1` and `KCNE3` remain more frequently detected in intestinal-like cells in all three groups. The effect becomes smaller after the ambiguous cells are removed, but its direction remains stable across patients and continuous state scores.

`KCNE2` follows the same detection direction, but the result is weaker and is not stable across all controls.

The full numbers and interpretation are in [RESULTS_AND_INTERPRETATION.md](RESULTS_AND_INTERPRETATION.md).

## Folder structure

```text
sensitivity_strict_malignancy/
|-- scripts/
|   `-- run_sensitivity_analysis.R
|-- outputs/
|   |-- figures/
|   `-- tables/
|-- README.md
`-- RESULTS_AND_INTERPRETATION.md
```

## Run

From this folder:

```powershell
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" scripts\run_sensitivity_analysis.R
```

The script reads `../outputs/objects/final_intestinal_diffuse_emt_seurat.rds`. It writes only inside this sensitivity folder.

Main figure: `outputs/figures/GSE183904_KCN_malignancy_sensitivity.pdf`

Main table: `outputs/tables/KCN_direction_stability.tsv`
