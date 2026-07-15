suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(Matrix)
  library(data.table)
  library(openxlsx)
  library(ggplot2)
  library(patchwork)
})

set.seed(12345)
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
gastric_root <- normalizePath(file.path(root, ".."), winslash = "/")
table_dir <- file.path(root, "outputs", "tables")
object_dir <- file.path(root, "outputs", "objects")
figure_dir <- file.path(root, "outputs", "intermediate", "figures")

zscore <- function(x) {
  value <- as.numeric(scale(x))
  value[!is.finite(value)] <- 0
  value
}
theme_tsne <- function() {
  theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11.5, hjust = 0.5),
      plot.subtitle = element_text(size = 8.5, hjust = 0.5),
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 7.5)
    )
}
single_plot <- function(plot) if (inherits(plot, "patchwork")) plot[[1]] else plot

epi <- readRDS(file.path(object_dir, "epithelial_candidates_seurat.rds"))
DefaultAssay(epi) <- "RNA"
cnv <- fread(file.path(table_dir, "infercnv_cell_scores.tsv"))
zhou_path <- file.path(gastric_root, "ZhouEtAlmarkergenes.xlsx")
stopifnot(file.exists(zhou_path))

zhou <- as.data.table(read.xlsx(zhou_path, sheet = "Table S4", startRow = 2))
zhou <- zhou[!is.na(Gene) & nzchar(Gene)]
zhou[, Gene := toupper(Gene)]
zhou_malignant <- unique(zhou[Cluster == "Malignant" & avg_logFC >= 0.50 & p_val_adj <= 0.01 & pct.1 >= 0.25, Gene])
zhou_nonmalignant <- unique(zhou[Cluster == "Nonmalignant" & avg_logFC >= 0.50 & p_val_adj <= 0.01 & pct.1 >= 0.25, Gene])

intestinal_genes <- c(
  "CDH17", "REG4", "MUC13", "TFF3", "CDX1", "CDX2", "KRT20", "FABP1",
  "LGALS4", "PIGR", "PHGR1", "GPA33", "MUC3A", "ANPEP", "AQP3", "FCGBP",
  "AGR2", "CES2", "ACSL5", "DGAT1", "TSPAN13", "SMIM24"
)
emt_up_genes <- c(
  "VIM", "CDH2", "SNAI2", "TWIST2", "ZEB1", "ZEB2", "FN1", "WNT5A",
  "NOTCH2", "AXL", "CALD1", "TPM1", "TPM2", "SERPINE2", "NNMT", "NUPR1",
  "TUSC3", "PRKD1", "GAS6", "THBS2", "CDH11", "LOXL1", "RECK", "DDR2",
  "EFEMP1", "CYP1B1", "LTBP2", "ROR2", "PDGFRB", "EDNRA", "ROBO1",
  "MEF2C", "IGFBP5", "S100A4", "TAGLN", "EGR1"
)
epithelial_junction_genes <- c(
  "EPCAM", "CDH1", "CLDN3", "CLDN4", "CLDN7", "KRT8", "KRT18", "KRT19",
  "MUC1", "TACSTD2"
)
signatures <- list(
  Zhou_malignant = zhou_malignant,
  Zhou_nonmalignant = zhou_nonmalignant,
  Intestinal = intestinal_genes,
  EMT_up = emt_up_genes,
  Epithelial_junction = epithelial_junction_genes
)
available <- lapply(signatures, intersect, y = rownames(epi))
if (any(lengths(available) < 4L)) {
  stop("Insufficient signature coverage: ", paste(names(available)[lengths(available) < 4L], collapse = ", "))
}
registry <- rbindlist(lapply(names(signatures), function(signature) {
  data.table(
    signature = signature,
    gene = signatures[[signature]],
    available = signatures[[signature]] %chin% rownames(epi),
    source = fcase(
      signature %chin% c("Zhou_malignant", "Zhou_nonmalignant"),
      "Zhou et al. 2023, Cellular and Molecular Life Sciences, Supplementary Table S4, DOI: 10.1007/s00018-023-04702-1",
      signature == "Intestinal",
      "Kim et al. 2022, npj Precision Oncology, GSE150290 intestinal tumor markers and gastric intestinal-lineage genes, DOI: 10.1038/s41698-022-00251-1",
      signature == "EMT_up",
      "Tanabe et al. 2014, International Journal of Oncology diffuse-type gastric cancer EMT signature, DOI: 10.3892/ijo.2014.2387; supported by Kim et al. 2022 GSE150290 diffuse/EMT controls",
      signature == "Epithelial_junction",
      "Canonical epithelial and junction markers used as EMT counter-score; motivated by Tanabe et al. 2014 and Kim et al. 2022"
    )
  )
}))
registry[, classifier_kcn_overlap := gene %chin% c("KCNQ1", "KCNE2", "KCNE3")]
fwrite(registry, file.path(table_dir, "final_signature_gene_registry.tsv"), sep = "\t")
stopifnot(!any(registry$classifier_kcn_overlap))

epi <- AddModuleScore_UCell(
  epi,
  features = available[c("Zhou_malignant", "Zhou_nonmalignant")],
  assay = "RNA",
  slot = "counts",
  maxRank = 1500,
  ncores = 1,
  missing_genes = "skip"
)
epi_meta <- as.data.table(epi@meta.data, keep.rownames = "cell_id")
epi_meta <- merge(epi_meta, cnv[, .(cell_id, cnv_score, cnv_high)], by = "cell_id", all.x = TRUE, sort = FALSE)
epi_meta[, malignant_program_z := zscore(Zhou_malignant_UCell)]
epi_meta[, nonmalignant_program_z := zscore(Zhou_nonmalignant_UCell)]
epi_meta[, tumor_program_delta := malignant_program_z - nonmalignant_program_z]
marker_threshold <- quantile(epi_meta[tissue == "Primary_Normal", tumor_program_delta], 0.95, na.rm = TRUE)
epi_meta[, marker_supported_tumor := tissue == "Primary_Tumor" & tumor_program_delta > marker_threshold]
epi_meta[, cnv_rescued_tumor := tissue == "Primary_Tumor" & !marker_supported_tumor & cnv_high %in% TRUE]
epi_meta[, malignant_final := marker_supported_tumor | cnv_rescued_tumor]
epi_meta[, selection_class := fifelse(
  marker_supported_tumor, "Tumor program",
  fifelse(cnv_rescued_tumor, "Added by inferCNV", "Not selected")
)]

epi_add <- as.data.frame(epi_meta[match(colnames(epi), cell_id), .(
  malignant_program_z, nonmalignant_program_z, tumor_program_delta, cnv_score, cnv_high,
  marker_supported_tumor, cnv_rescued_tumor, malignant_final, selection_class
)])
rownames(epi_add) <- colnames(epi)
epi <- AddMetaData(epi, epi_add)
saveRDS(epi, file.path(object_dir, "epithelial_with_malignancy_scores_seurat.rds"), compress = FALSE)
fwrite(epi_meta, file.path(table_dir, "epithelial_cells_with_malignancy_scores.tsv.gz"), sep = "\t", compress = "gzip")

selection_summary <- epi_meta[, .(
  cells = .N,
  median_tumor_delta = median(tumor_program_delta),
  median_cnv_score = median(cnv_score, na.rm = TRUE)
), by = .(tissue, selection_class)][order(tissue, selection_class)]
fwrite(selection_summary, file.path(table_dir, "malignant_selection_summary.tsv"), sep = "\t")
fwrite(epi_meta[tissue == "Primary_Tumor", .N, by = .(patient, lauren, selection_class)],
       file.path(table_dir, "malignant_selection_by_patient.tsv"), sep = "\t")
fwrite(data.table(parameter = "normal_tumor_program_delta_q95", value = marker_threshold),
       file.path(table_dir, "malignancy_thresholds.tsv"), sep = "\t")

malignant_cells <- epi_meta[malignant_final == TRUE, cell_id]
tumor_pre_qc <- subset(epi, cells = malignant_cells)
immune_markers <- intersect(
  c("PTPRC", "CD3D", "CD3E", "NKG7", "LST1", "TYROBP", "FCER1G", "CD52", "LCP2", "SRGN"),
  rownames(tumor_pre_qc)
)
pre_qc_counts <- GetAssayData(tumor_pre_qc, assay = "RNA", layer = "counts")
immune_detected <- Matrix::colSums(pre_qc_counts[immune_markers, , drop = FALSE] > 0)
contamination <- data.table(
  cell_id = names(immune_detected),
  immune_markers_detected = as.integer(immune_detected),
  removed_as_immune_contaminant = immune_detected >= 2L
)
contamination <- merge(contamination, epi_meta[, .(cell_id, patient, lauren, selection_class)],
                       by = "cell_id", all.x = TRUE, sort = FALSE)
fwrite(contamination, file.path(table_dir, "immune_contamination_filter.tsv"), sep = "\t")
keep_cells <- contamination[removed_as_immune_contaminant == FALSE, cell_id]
tumor <- subset(tumor_pre_qc, cells = keep_cells)

tumor <- NormalizeData(tumor, verbose = FALSE)
tumor <- FindVariableFeatures(tumor, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
tumor <- ScaleData(tumor, vars.to.regress = "nCount_RNA", features = VariableFeatures(tumor), verbose = FALSE)
tumor <- RunPCA(tumor, features = VariableFeatures(tumor), npcs = 30, seed.use = 12345, verbose = FALSE)
tumor <- FindNeighbors(tumor, dims = 1:20, k.param = 30, verbose = FALSE)
tumor <- FindClusters(tumor, resolution = 0.8, random.seed = 12345, verbose = FALSE)
tumor <- RunTSNE(tumor, dims = 1:20, seed.use = 12345, check_duplicates = FALSE, verbose = FALSE)
tumor <- AddModuleScore_UCell(
  tumor,
  features = available[c("Intestinal", "EMT_up", "Epithelial_junction")],
  assay = "RNA",
  slot = "counts",
  maxRank = 1500,
  ncores = 1,
  missing_genes = "skip"
)

tumor_meta <- as.data.table(tumor@meta.data, keep.rownames = "cell_id")
tumor_meta[, intestinal_z := zscore(Intestinal_UCell)]
tumor_meta[, emt_up_z := zscore(EMT_up_UCell)]
tumor_meta[, epithelial_junction_z := zscore(Epithelial_junction_UCell)]
tumor_meta[, diffuse_emt_raw := emt_up_z - epithelial_junction_z]
tumor_meta[, diffuse_emt_z := zscore(diffuse_emt_raw)]
tumor_meta[, state_delta := intestinal_z - diffuse_emt_z]
confidence_margin <- 0.50
tumor_meta[, state_final := fifelse(
  state_delta >= confidence_margin & intestinal_z >= 0,
  "Intestinal-like",
  fifelse(state_delta <= -confidence_margin & diffuse_emt_z >= 0, "Diffuse/EMT-like", "Indeterminate")
)]
tumor_add <- as.data.frame(tumor_meta[match(colnames(tumor), cell_id), .(
  intestinal_z, emt_up_z, epithelial_junction_z, diffuse_emt_z, state_delta, state_final
)])
rownames(tumor_add) <- colnames(tumor)
tumor <- AddMetaData(tumor, tumor_add)
saveRDS(tumor, file.path(object_dir, "malignant_unbiased_states_seurat.rds"), compress = FALSE)
fwrite(tumor_meta, file.path(table_dir, "malignant_cells_with_states.tsv.gz"), sep = "\t", compress = "gzip")

final_cells <- tumor_meta[state_final %chin% c("Intestinal-like", "Diffuse/EMT-like"), cell_id]
final_obj <- subset(tumor, cells = final_cells)
saveRDS(final_obj, file.path(object_dir, "final_intestinal_diffuse_emt_seurat.rds"), compress = FALSE)

state_summary <- tumor_meta[, .(
  cells = .N,
  patients = uniqueN(patient),
  median_intestinal_z = median(intestinal_z),
  median_diffuse_emt_z = median(diffuse_emt_z),
  median_state_delta = median(state_delta)
), by = state_final][order(state_final)]
fwrite(state_summary, file.path(table_dir, "state_summary.tsv"), sep = "\t")
state_by_patient <- tumor_meta[, .N, by = .(patient, lauren, state_final)]
state_by_patient[, patient_total := sum(N), by = patient]
state_by_patient[, fraction := N / patient_total]
fwrite(state_by_patient, file.path(table_dir, "state_composition_by_patient.tsv"), sep = "\t")
sensitivity <- rbindlist(lapply(c(0, 0.25, 0.5, 0.75, 1), function(margin) {
  tumor_meta[, .(
    margin = margin,
    retained_cells = sum((state_delta >= margin & intestinal_z >= 0) |
                         (state_delta <= -margin & diffuse_emt_z >= 0)),
    intestinal_cells = sum(state_delta >= margin & intestinal_z >= 0),
    diffuse_emt_cells = sum(state_delta <= -margin & diffuse_emt_z >= 0)
  )]
}))
sensitivity[, retained_fraction := retained_cells / nrow(tumor_meta)]
fwrite(sensitivity, file.path(table_dir, "state_confidence_sensitivity.tsv"), sep = "\t")

p1 <- FeaturePlot(epi, "tumor_program_delta", reduction = "tsne", pt.size = 0.16,
                  order = TRUE, cols = c("#ffd23f", "#ff6b6b")) +
  labs(title = "Malignant epithelial program", subtitle = "Published malignant minus non-malignant score") + theme_tsne()
p2 <- FeaturePlot(epi, "cnv_score", reduction = "tsne", pt.size = 0.16,
                  order = TRUE, cols = c("#88d498", "#74b9ff")) +
  labs(title = "inferCNV support", subtitle = "Normal epithelial reference") + theme_tsne()
p3 <- DimPlot(epi, reduction = "tsne", group.by = "selection_class", pt.size = 0.16,
              order = c("Not selected", "Added by inferCNV", "Tumor program"),
              cols = c("Not selected" = "#b8a9fa", "Added by inferCNV" = "#74b9ff", "Tumor program" = "#ffa552")) +
  labs(title = "Malignant-cell selection") + theme_tsne()
p4 <- DimPlot(tumor, reduction = "tsne", group.by = "patient", pt.size = 0.18) +
  labs(title = "Cleaned malignant cells", subtitle = "Tumor-only t-SNE by patient") + theme_tsne()
page_malignancy <- wrap_plots(lapply(list(p1, p2, p3, p4), single_plot), ncol = 2) +
  plot_annotation(title = "Malignant epithelial-cell selection")

p5 <- FeaturePlot(tumor, "intestinal_z", reduction = "tsne", pt.size = 0.20,
                  order = TRUE, cols = c("#f2f2f2", "#2b9348")) +
  labs(title = "Intestinal program", subtitle = "UCell z-score") + theme_tsne()
p6 <- FeaturePlot(tumor, "diffuse_emt_z", reduction = "tsne", pt.size = 0.20,
                  order = TRUE, cols = c("#f2f2f2", "#ff5714")) +
  labs(title = "Diffuse/EMT program", subtitle = "EMT activation minus epithelial-junction preservation") + theme_tsne()
p7 <- DimPlot(tumor, reduction = "tsne", group.by = "state_final", pt.size = 0.20,
              order = c("Indeterminate", "Diffuse/EMT-like", "Intestinal-like"),
              cols = c("Indeterminate" = "#b8a9fa", "Diffuse/EMT-like" = "#ff5714", "Intestinal-like" = "#2b9348")) +
  labs(title = "High-confidence state assignment") + theme_tsne()
p8 <- DimPlot(final_obj, reduction = "tsne", group.by = "state_final", pt.size = 0.22,
              cols = c("Diffuse/EMT-like" = "#ff5714", "Intestinal-like" = "#2b9348")) +
  labs(title = "Final two-state object", subtitle = paste0(ncol(final_obj), " cells")) + theme_tsne()
page_states <- wrap_plots(lapply(list(p5, p6, p7, p8), single_plot), ncol = 2) +
  plot_annotation(title = "Intestinal and diffuse/EMT malignant-cell states")

pdf(file.path(figure_dir, "02_malignancy_and_states_tsne.pdf"), width = 10, height = 8,
    onefile = TRUE, useDingbats = FALSE)
print(page_malignancy)
print(page_states)
dev.off()

initial_flow <- fread(file.path(table_dir, "cell_flow_initial.tsv"))
flow <- rbind(
  initial_flow,
  data.table(step = c("Malignant candidates", "After immune-contamination filter", "Final high-confidence states"),
             cells = c(length(malignant_cells), ncol(tumor), ncol(final_obj)))
)
fwrite(flow, file.path(table_dir, "cell_flow_complete.tsv"), sep = "\t")
message("State pipeline completed: ", ncol(epi), " epithelial candidates; ",
        length(malignant_cells), " malignant candidates; ", ncol(tumor),
        " cleaned malignant cells; ", ncol(final_obj), " final cells.")
