suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

set.seed(12345)
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
analysis_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
dataset_root <- normalizePath(file.path(analysis_root, ".."), winslash = "/")
table_dir <- file.path(analysis_root, "outputs", "tables")
figure_dir <- file.path(analysis_root, "outputs", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

object_path <- file.path(dataset_root, "outputs", "objects", "final_intestinal_diffuse_emt_seurat.rds")
obj <- readRDS(object_path)
DefaultAssay(obj) <- "RNA"
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
expression <- GetAssayData(obj, assay = "RNA", layer = "data")
meta <- as.data.table(obj@meta.data, keep.rownames = "cell_id")

kcn_genes <- c("KCNQ1", "KCNE2", "KCNE3")
stopifnot(all(kcn_genes %in% rownames(counts)))

epithelial_core <- intersect(c("KRT8", "KRT18", "KRT19", "EPCAM", "CDH1", "MUC1", "TACSTD2"), rownames(counts))
fibroblast_core <- intersect(c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FN1", "SPARC"), rownames(counts))
mesothelial_core <- intersect(c("CALB2", "WT1", "UPK3B", "LRRN4", "MSLN"), rownames(counts))

detected_count <- function(genes) Matrix::colSums(counts[genes, meta$cell_id, drop = FALSE] > 0)
meta[, epithelial_core_detected := as.integer(detected_count(epithelial_core))]
meta[, fibroblast_core_detected := as.integer(detected_count(fibroblast_core))]
meta[, mesothelial_core_detected := as.integer(detected_count(mesothelial_core))]

ambiguous_clusters <- c("2", "5", "16", "17", "18")
meta[, original_final := TRUE]
meta[, cluster_screened := !as.character(seurat_clusters) %chin% ambiguous_clusters]
meta[, orthogonal_strict :=
  tissue == "Primary_Tumor" &
  epithelial_core_detected >= 2L &
  mesothelial_core_detected < 2L &
  (
    cnv_high %in% TRUE |
    (
      marker_supported_tumor %in% TRUE &
      epithelial_score > fibroblast_score &
      fibroblast_core_detected < 3L &
      epithelial_core_detected >= 4L
    )
  )
]

cohort_columns <- c(
  "Original final" = "original_final",
  "Cluster-screened" = "cluster_screened",
  "Orthogonal strict" = "orthogonal_strict"
)

parameters <- data.table(
  parameter = c(
    "input_object", "ambiguous_clusters", "epithelial_core", "fibroblast_core", "mesothelial_core",
    "strict_min_epithelial_genes", "strict_max_mesothelial_genes", "marker_branch_min_epithelial_genes",
    "marker_branch_max_fibroblast_genes", "minimum_cells_per_state_for_patient_comparison"
  ),
  value = c(
    "../outputs/objects/final_intestinal_diffuse_emt_seurat.rds",
    paste(ambiguous_clusters, collapse = ","), paste(epithelial_core, collapse = ","),
    paste(fibroblast_core, collapse = ","), paste(mesothelial_core, collapse = ","),
    "2", "1", "4", "2", "10"
  )
)
fwrite(parameters, file.path(table_dir, "analysis_parameters.tsv"), sep = "\t")

cohort_membership <- rbindlist(lapply(names(cohort_columns), function(cohort) {
  keep <- meta[[cohort_columns[[cohort]]]]
  meta[keep, .(cell_id, patient, state = state_final, seurat_clusters)][, cohort := cohort]
}))
setcolorder(cohort_membership, c("cohort", "cell_id", "patient", "state", "seurat_clusters"))
fwrite(cohort_membership, file.path(table_dir, "cohort_cell_membership.tsv.gz"), sep = "\t", compress = "gzip")

cohort_summary <- cohort_membership[, .(
  cells = .N,
  patients = uniqueN(patient)
), by = .(cohort, state)]
cohort_summary[, cohort_total := sum(cells), by = cohort]
cohort_summary[, state_fraction := cells / cohort_total]
fwrite(cohort_summary, file.path(table_dir, "cohort_summary.tsv"), sep = "\t")

cluster_audit <- meta[, .(
  cells = .N,
  patients = uniqueN(patient),
  diffuse_emt_fraction = mean(state_final == "Diffuse/EMT-like"),
  cnv_high_fraction = mean(cnv_high %in% TRUE),
  marker_supported_fraction = mean(marker_supported_tumor %in% TRUE),
  median_epithelial_core = as.numeric(median(epithelial_core_detected)),
  median_fibroblast_core = as.numeric(median(fibroblast_core_detected)),
  median_mesothelial_core = as.numeric(median(mesothelial_core_detected)),
  epithelial_dominant_fraction = mean(epithelial_score > fibroblast_score),
  excluded_cluster_screen = as.character(first(seurat_clusters)) %chin% ambiguous_clusters
), by = seurat_clusters][order(as.integer(as.character(seurat_clusters)))]
fwrite(cluster_audit, file.path(table_dir, "cluster_orthogonal_audit.tsv"), sep = "\t")

detection_cell <- rbindlist(lapply(names(cohort_columns), function(cohort) {
  ids <- meta[get(cohort_columns[[cohort]]), cell_id]
  cell_meta <- meta[match(ids, cell_id)]
  rbindlist(lapply(kcn_genes, function(gene) {
    data.table(
      cohort = cohort,
      cell_id = ids,
      patient = cell_meta$patient,
      state = cell_meta$state_final,
      gene = gene,
      raw_count = as.numeric(counts[gene, ids])
    )[, detected := raw_count > 0]
  }))
}))

detection <- detection_cell[, .(
  cells = .N,
  expressing_cells = sum(detected),
  detection_fraction = mean(detected),
  raw_transcripts = sum(raw_count),
  mean_raw_count_all_cells = mean(raw_count),
  mean_raw_count_expressing_cells = if (any(detected)) mean(raw_count[detected]) else NA_real_
), by = .(cohort, gene, state)]
fwrite(detection, file.path(table_dir, "KCN_detection_by_cohort_state.tsv"), sep = "\t")

fisher_results <- rbindlist(lapply(names(cohort_columns), function(cohort_id) {
  rbindlist(lapply(kcn_genes, function(gene_id) {
    rows <- detection[cohort == cohort_id & gene == gene_id]
    tab <- matrix(c(
      rows[state == "Intestinal-like", expressing_cells],
      rows[state == "Intestinal-like", cells - expressing_cells],
      rows[state == "Diffuse/EMT-like", expressing_cells],
      rows[state == "Diffuse/EMT-like", cells - expressing_cells]
    ), nrow = 2, byrow = TRUE)
    test <- fisher.test(tab)
    data.table(cohort = cohort_id, gene = gene_id, odds_ratio = unname(test$estimate), p_value = test$p.value)
  }))
}))
fisher_results[, FDR := p.adjust(p_value, method = "BH"), by = cohort]
fisher_results <- merge(
  fisher_results,
  dcast(detection, cohort + gene ~ state, value.var = "detection_fraction"),
  by = c("cohort", "gene")
)
fisher_results[, detection_difference := `Intestinal-like` - `Diffuse/EMT-like`]
fwrite(fisher_results, file.path(table_dir, "KCN_fisher_sensitivity.tsv"), sep = "\t")

patient_detection <- detection_cell[, .(
  cells = .N,
  detection_fraction = mean(detected)
), by = .(cohort, gene, patient, state)]
paired <- dcast(patient_detection, cohort + gene + patient ~ state, value.var = c("cells", "detection_fraction"))
paired <- paired[
  !is.na(`detection_fraction_Intestinal-like`) & !is.na(`detection_fraction_Diffuse/EMT-like`) &
    `cells_Intestinal-like` >= 10L & `cells_Diffuse/EMT-like` >= 10L
]
paired[, detection_difference := `detection_fraction_Intestinal-like` - `detection_fraction_Diffuse/EMT-like`]
patient_direction <- paired[, {
  test <- if (.N >= 3L && any(detection_difference != 0)) wilcox.test(detection_difference, mu = 0, exact = FALSE) else NULL
  .(
    eligible_patients = .N,
    median_patient_difference = median(detection_difference),
    patients_positive = sum(detection_difference > 0),
    patients_negative = sum(detection_difference < 0),
    wilcoxon_p_value = if (is.null(test)) NA_real_ else test$p.value
  )
}, by = .(cohort, gene)]
patient_direction[, FDR := p.adjust(wilcoxon_p_value, method = "BH"), by = cohort]
fwrite(paired, file.path(table_dir, "KCN_patient_detection_differences.tsv"), sep = "\t")
fwrite(patient_direction, file.path(table_dir, "KCN_patient_direction_summary.tsv"), sep = "\t")

within_patient_rho <- function(x, y, patient) {
  keep <- is.finite(x) & is.finite(y) & !is.na(patient)
  x <- x[keep]
  y <- y[keep]
  patient <- patient[keep]
  eligible <- names(which(table(patient) >= 20L))
  keep <- patient %in% eligible
  x <- x[keep]
  y <- y[keep]
  patient <- patient[keep]
  groups <- split(seq_along(x), patient)
  xr <- yr <- numeric(length(x))
  for (index in groups) {
    xr[index] <- rank(x[index], ties.method = "average")
    yr[index] <- rank(y[index], ties.method = "average")
    xr[index] <- xr[index] - mean(xr[index])
    yr[index] <- yr[index] - mean(yr[index])
  }
  data.table(rho = cor(xr, yr), cells = length(x), patients = length(groups))
}

score_columns <- c(Intestinal = "intestinal_z", `Diffuse/EMT` = "diffuse_emt_z")
correlations <- rbindlist(lapply(names(cohort_columns), function(cohort) {
  ids <- meta[get(cohort_columns[[cohort]]), cell_id]
  cohort_meta <- meta[match(ids, cell_id)]
  rbindlist(lapply(kcn_genes, function(gene) {
    gene_expression <- as.numeric(expression[gene, ids])
    rbindlist(lapply(names(score_columns), function(signature) {
      score <- cohort_meta[[score_columns[[signature]]]]
      pooled <- suppressWarnings(cor.test(gene_expression, score, method = "spearman", exact = FALSE))
      within <- within_patient_rho(gene_expression, score, cohort_meta$patient)
      rbind(
        data.table(cohort = cohort, gene = gene, signature = signature, method = "Pooled Spearman",
                   rho = unname(pooled$estimate), p_value = pooled$p.value,
                   cells = length(ids), patients = uniqueN(cohort_meta$patient)),
        data.table(cohort = cohort, gene = gene, signature = signature, method = "Within-patient ranks",
                   rho = within$rho, p_value = NA_real_, cells = within$cells, patients = within$patients)
      )
    }))
  }))
}))
correlations[, FDR := p.adjust(p_value, method = "BH"), by = .(cohort, method)]
fwrite(correlations, file.path(table_dir, "KCN_signature_correlation_sensitivity.tsv"), sep = "\t")

direction_stability <- fisher_results[, .(
  cohorts = .N,
  all_detection_differences_intestinal_positive = all(detection_difference > 0),
  minimum_detection_difference = min(detection_difference),
  maximum_detection_difference = max(detection_difference),
  minimum_odds_ratio = min(odds_ratio),
  maximum_FDR = max(FDR)
), by = gene]
within_stability <- correlations[method == "Within-patient ranks", .(
  intestinal_rho_min = min(rho[signature == "Intestinal"]),
  intestinal_rho_max = max(rho[signature == "Intestinal"]),
  diffuse_rho_min = min(rho[signature == "Diffuse/EMT"]),
  diffuse_rho_max = max(rho[signature == "Diffuse/EMT"])
), by = gene]
direction_stability <- merge(direction_stability, within_stability, by = "gene")
direction_stability[, concordant_signature_direction := intestinal_rho_min > 0 & diffuse_rho_max < 0]
fwrite(direction_stability, file.path(table_dir, "KCN_direction_stability.tsv"), sep = "\t")

cohort_levels <- names(cohort_columns)
detection[, cohort := factor(cohort, levels = cohort_levels)]
fisher_results[, cohort := factor(cohort, levels = cohort_levels)]
correlations[, cohort := factor(cohort, levels = cohort_levels)]
state_colors <- c("Intestinal-like" = "#2b9348", "Diffuse/EMT-like" = "#ff5714")

p_counts <- ggplot(cohort_summary, aes(cohort, cells, fill = state)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = cells), position = position_stack(vjust = 0.5), size = 3.3, fontface = "bold") +
  scale_fill_manual(values = state_colors) +
  labs(title = "Cells retained by sensitivity cohort", x = NULL, y = "Cells", fill = "State") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 18, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "top")

p_detection <- ggplot(detection, aes(cohort, 100 * detection_fraction, fill = state)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.64) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * detection_fraction)),
            position = position_dodge(width = 0.75), vjust = -0.25, size = 2.8) +
  facet_wrap(~gene, nrow = 1) +
  scale_fill_manual(values = state_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "KCN detection remains directional after stricter controls", x = NULL,
       y = "Expressing cells (%)", fill = "State") +
  theme_bw(base_size = 10.5) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 18, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "top")

correlations[, method := factor(method, levels = c("Pooled Spearman", "Within-patient ranks"))]
p_cor <- ggplot(correlations, aes(rho, cohort, color = method, shape = method)) +
  geom_vline(xintercept = 0, color = "grey70") +
  geom_point(size = 2.7, position = position_dodge(width = 0.35)) +
  facet_grid(gene ~ signature) +
  scale_color_manual(values = c("Pooled Spearman" = "#74b9ff", "Within-patient ranks" = "#6247aa")) +
  labs(title = "Continuous signature associations", x = expression(rho), y = NULL, color = NULL, shape = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "top", strip.background = element_blank(), strip.text = element_text(face = "bold"))

pdf(file.path(figure_dir, "GSE183904_KCN_malignancy_sensitivity.pdf"), width = 10, height = 8, onefile = TRUE, useDingbats = FALSE)
print(p_counts)
print(p_detection)
print(p_cor)
dev.off()

message("Sensitivity analysis completed without modifying the main pipeline.")
