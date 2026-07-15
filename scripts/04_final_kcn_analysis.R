suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(qpdf)
})

set.seed(12345)
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
object_dir <- file.path(root, "outputs", "objects")
source_figure_dir <- file.path(root, "outputs", "intermediate", "figures")
final_table_dir <- file.path(root, "outputs", "final", "tables")
final_figure_dir <- file.path(root, "outputs", "final", "figures")
dir.create(final_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_figure_dir, recursive = TRUE, showWarnings = FALSE)

intestinal_color <- "#2b9348"
diffuse_color <- "#ff5714"
kcn_expressed <- "#6247aa"
kcn_other <- "#ffd500"
state_colors <- c("Intestinal-like" = intestinal_color, "Diffuse/EMT-like" = diffuse_color)

theme_tsne <- function() {
  theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11.5, hjust = 0.5),
      plot.subtitle = element_text(size = 8.5, hjust = 0.5),
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8)
    )
}
single_plot <- function(plot) if (inherits(plot, "patchwork")) plot[[1]] else plot

obj <- readRDS(file.path(object_dir, "final_intestinal_diffuse_emt_seurat.rds"))
DefaultAssay(obj) <- "RNA"
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
expression <- GetAssayData(obj, assay = "RNA", layer = "data")
kcn_genes <- intersect(c("KCNQ1", "KCNE2", "KCNE3"), rownames(counts))
stopifnot(length(kcn_genes) == 3L)
meta <- as.data.table(obj@meta.data, keep.rownames = "cell_id")

for (gene in kcn_genes) {
  status <- ifelse(as.numeric(counts[gene, colnames(obj)]) > 0, "Expressed", "Not expressed")
  obj[[paste0(gene, "_status")]] <- factor(status, levels = c("Not expressed", "Expressed"))
}

detection_cell <- rbindlist(lapply(kcn_genes, function(gene) {
  data.table(
    cell_id = meta$cell_id,
    patient = meta$patient,
    state = meta$state_final,
    gene = gene,
    raw_count = as.numeric(counts[gene, meta$cell_id])
  )[, detected := raw_count > 0]
}))
detection <- detection_cell[, .(
  cells = .N,
  expressing_cells = sum(detected),
  non_expressing_cells = sum(!detected),
  detection_fraction = mean(detected),
  raw_transcripts = sum(raw_count),
  mean_raw_count_all_cells = mean(raw_count),
  mean_raw_count_expressing_cells = ifelse(any(detected), mean(raw_count[detected]), NA_real_)
), by = .(gene, state)]
patient_detection <- detection_cell[, .(
  cells = .N,
  expressing_cells = sum(detected),
  detection_fraction = mean(detected),
  raw_transcripts = sum(raw_count)
), by = .(gene, patient, state)]

fisher_results <- rbindlist(lapply(kcn_genes, function(gene_id) {
  rows <- detection[gene == gene_id]
  matrix_2x2 <- matrix(c(
    rows[state == "Intestinal-like", expressing_cells], rows[state == "Intestinal-like", non_expressing_cells],
    rows[state == "Diffuse/EMT-like", expressing_cells], rows[state == "Diffuse/EMT-like", non_expressing_cells]
  ), nrow = 2, byrow = TRUE,
  dimnames = list(state = c("Intestinal-like", "Diffuse/EMT-like"), detection = c("Expressed", "Not expressed")))
  test <- fisher.test(matrix_2x2)
  data.table(
    gene = gene_id,
    test = "Fisher exact test",
    patient_stratified = FALSE,
    patients = uniqueN(meta$patient),
    odds_ratio = unname(test$estimate),
    p_value = test$p.value
  )
}))

detection_tests <- fisher_results
detection_tests[, FDR := p.adjust(p_value, method = "BH")]
detection_tests <- merge(
  detection_tests,
  dcast(detection, gene ~ state, value.var = "detection_fraction"),
  by = "gene", all.x = TRUE
)
detection_tests[, detection_difference := `Intestinal-like` - `Diffuse/EMT-like`]

within_patient_rank_test <- function(x, y, patient, permutations = 1000L) {
  keep <- is.finite(x) & is.finite(y) & !is.na(patient)
  x <- x[keep]
  y <- y[keep]
  patient <- patient[keep]
  eligible <- names(which(table(patient) >= 20L))
  keep <- patient %in% eligible
  x <- x[keep]
  y <- y[keep]
  patient <- patient[keep]
  split_index <- split(seq_along(x), patient)
  x_rank <- y_rank <- numeric(length(x))
  for (index in split_index) {
    x_rank[index] <- rank(x[index], ties.method = "average")
    y_rank[index] <- rank(y[index], ties.method = "average")
    x_rank[index] <- x_rank[index] - mean(x_rank[index])
    y_rank[index] <- y_rank[index] - mean(y_rank[index])
  }
  observed <- cor(x_rank, y_rank)
  null <- numeric(permutations)
  for (iteration in seq_len(permutations)) {
    permuted <- x_rank
    for (index in split_index) permuted[index] <- sample(x_rank[index])
    null[iteration] <- cor(permuted, y_rank)
  }
  list(rho = observed, p_value = (1 + sum(abs(null) >= abs(observed))) / (permutations + 1),
       cells = length(x), patients = length(split_index))
}

score_columns <- c(Intestinal = "intestinal_z", `Diffuse/EMT` = "diffuse_emt_z")
correlations <- rbindlist(lapply(kcn_genes, function(gene) {
  gene_expression <- as.numeric(expression[gene, meta$cell_id])
  rbindlist(lapply(names(score_columns), function(signature) {
    score <- meta[[score_columns[[signature]]]]
    pooled <- suppressWarnings(cor.test(gene_expression, score, method = "spearman", exact = FALSE))
    within <- within_patient_rank_test(gene_expression, score, meta$patient, permutations = 1000L)
    rbind(
      data.table(gene = gene, signature = signature, method = "Pooled Spearman",
                 patient_stratified = FALSE, cells = length(gene_expression), patients = uniqueN(meta$patient),
                 rho = unname(pooled$estimate), p_value = pooled$p.value),
      data.table(gene = gene, signature = signature, method = "Within-patient rank permutation",
                 patient_stratified = TRUE, cells = within$cells, patients = within$patients,
                 rho = within$rho, p_value = within$p_value)
    )
  }))
}))
correlations[, FDR := p.adjust(p_value, method = "BH"), by = method]

fwrite(detection, file.path(final_table_dir, "KCN_detection_by_state.tsv"), sep = "\t")
fwrite(patient_detection, file.path(final_table_dir, "KCN_detection_by_patient_state.tsv"), sep = "\t")
fwrite(detection_tests, file.path(final_table_dir, "KCN_detection_tests.tsv"), sep = "\t")
fwrite(correlations, file.path(final_table_dir, "KCN_signature_correlations.tsv"), sep = "\t")
fwrite(meta, file.path(final_table_dir, "final_cell_metadata.tsv.gz"), sep = "\t", compress = "gzip")

patient_state <- meta[, .N, by = .(patient, lauren, state_final)]
patient_state_wide <- dcast(patient_state, patient + lauren ~ state_final, value.var = "N", fill = 0)
patient_state_wide[, final_cells := `Intestinal-like` + `Diffuse/EMT-like`]
patient_state_wide[, intestinal_fraction := `Intestinal-like` / final_cells]
patient_state_wide[, diffuse_emt_fraction := `Diffuse/EMT-like` / final_cells]
fwrite(patient_state_wide, file.path(final_table_dir, "final_state_composition_by_patient.tsv"), sep = "\t")

lauren_check <- wilcox.test(diffuse_emt_fraction ~ lauren, data = patient_state_wide, exact = FALSE)
lauren_validation <- data.table(
  validation = "Diffuse/EMT-like fraction by patient Lauren type",
  test = "Wilcoxon rank-sum test",
  diffuse_patients = patient_state_wide[lauren == "Diffuse", .N],
  intestinal_patients = patient_state_wide[lauren == "Intestinal", .N],
  median_fraction_diffuse_patients = patient_state_wide[lauren == "Diffuse", median(diffuse_emt_fraction)],
  median_fraction_intestinal_patients = patient_state_wide[lauren == "Intestinal", median(diffuse_emt_fraction)],
  p_value = lauren_check$p.value
)
fwrite(lauren_validation, file.path(final_table_dir, "state_classifier_lauren_validation.tsv"), sep = "\t")

summary_table <- merge(
  detection,
  detection_tests[, .(gene, fisher_odds_ratio = odds_ratio, fisher_p_value = p_value, fisher_FDR = FDR)],
  by = "gene", all.x = TRUE
)
fwrite(summary_table, file.path(final_table_dir, "KCN_final_summary.tsv"), sep = "\t")

p_state <- DimPlot(obj, reduction = "tsne", group.by = "state_final", pt.size = 0.20,
                   cols = state_colors) + labs(title = "Final malignant-cell states") + theme_tsne()
plot_detection <- function(gene) {
  DimPlot(obj, reduction = "tsne", group.by = paste0(gene, "_status"), pt.size = 0.20,
          order = c("Not expressed", "Expressed"),
          cols = c("Not expressed" = kcn_other, "Expressed" = kcn_expressed)) +
    labs(title = gene, subtitle = "Raw-count detection", color = NULL) + theme_tsne()
}
page_projection <- wrap_plots(
  lapply(c(list(p_state), lapply(kcn_genes, plot_detection)), single_plot), ncol = 2
) + plot_annotation(title = "KCN detection in final intestinal and diffuse/EMT cells")

fisher_labels <- detection_tests[, .(gene, label = sprintf("Fisher FDR %.2g", FDR))]
detection_plot_data <- merge(detection, fisher_labels, by = "gene", all.x = TRUE)
p_detection <- ggplot(detection_plot_data, aes(gene, detection_fraction * 100, fill = state)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.68) +
  geom_text(aes(label = sprintf("%.1f%%", detection_fraction * 100)),
            position = position_dodge(width = 0.78), vjust = -0.25, size = 3.3, fontface = "bold") +
  geom_text(data = unique(detection_plot_data[, .(gene, label)]),
            aes(gene, Inf, label = label), inherit.aes = FALSE, vjust = 1.4, size = 3.2) +
  scale_fill_manual(values = state_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Detection rate in final cells", subtitle = "Fisher exact test on expressing versus non-expressing cells",
       x = NULL, y = "Expressing cells (%)", fill = "State") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5), panel.grid.minor = element_blank(),
        legend.position = "top")

correlations[, method_label := fifelse(patient_stratified, "Within-patient", "Pooled")]
p_correlation <- ggplot(correlations, aes(rho, gene, color = method_label, shape = method_label)) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.5) +
  geom_point(size = 3, position = position_dodge(width = 0.45)) +
  facet_wrap(~ signature, nrow = 1) +
  scale_color_manual(values = c("Pooled" = "#74b9ff", "Within-patient" = "#6247aa")) +
  coord_cartesian(xlim = c(-0.5, 0.5)) +
  labs(title = "Association with continuous signatures", subtitle = "Spearman and within-patient rank analyses",
       x = expression(rho), y = NULL, color = NULL, shape = NULL) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5), panel.grid.minor = element_blank(),
        legend.position = "top", strip.background = element_blank(), strip.text = element_text(face = "bold"))

page_statistics <- p_detection / p_correlation +
  plot_annotation(title = "KCN detection and signature associations")
final_kcn_pdf <- file.path(final_figure_dir, "GSE183904_final_KCN_results.pdf")
pdf(final_kcn_pdf, width = 10, height = 8, onefile = TRUE, useDingbats = FALSE)
print(page_projection)
print(page_statistics)
dev.off()

complete_pdf <- file.path(final_figure_dir, "GSE183904_complete_pipeline_and_KCN_results.pdf")
qpdf::pdf_combine(
  input = c(
    file.path(source_figure_dir, "01_epithelial_preparation_tsne.pdf"),
    file.path(source_figure_dir, "02_malignancy_and_states_tsne.pdf"),
    final_kcn_pdf
  ),
  output = complete_pdf
)
message("Final KCN analysis completed for ", ncol(obj), " cells.")
