# Plot non-signature marker projections on the final dataset4 t-SNE.
# These genes were not used in the intestinal/diffuse-EMT state classifier.

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

root <- "C:/Users/pcyou/Desktop/stage_LLM/M1_Bioinformatics_Ion_Channel/Gastric_KCNQ1/gastric_dataset4"
rds_path <- file.path(root, "outputs/objects/final_intestinal_diffuse_emt_seurat.rds")
figure_dir <- file.path(root, "outputs/final/figures")
table_dir <- file.path(root, "outputs/final/tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(rds_path)
DefaultAssay(obj) <- "RNA"

registry <- fread(file.path(root, "outputs/tables/final_signature_gene_registry.tsv"))
state_signature_genes <- registry[
  signature %chin% c("Intestinal", "EMT_up", "Epithelial_junction"),
  unique(gene)
]

intestinal_markers <- c("VIL1", "CDHR5", "SPINK4", "MUC2", "APOA4")
diffuse_markers <- c("SPARC", "COL1A1", "COL3A1", "SERPINE1", "MMP2")
marker_table <- data.table(
  gene = c(intestinal_markers, diffuse_markers),
  marker_axis = c(rep("Intestinal independent", length(intestinal_markers)),
                  rep("Diffuse/EMT independent", length(diffuse_markers)))
)
marker_table[, available := gene %in% rownames(obj)]
marker_table[, used_in_state_signature := gene %in% state_signature_genes]

if (any(!marker_table$available)) {
  stop("Missing marker(s) in final object: ", paste(marker_table[available == FALSE, gene], collapse = ", "))
}
if (any(marker_table$used_in_state_signature)) {
  stop("Marker(s) were used in the state signatures: ",
       paste(marker_table[used_in_state_signature == TRUE, gene], collapse = ", "))
}

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
meta <- as.data.table(obj@meta.data, keep.rownames = "cell_id")

marker_metrics <- rbindlist(lapply(marker_table$gene, function(gene_id) {
  raw <- as.numeric(counts[gene_id, meta$cell_id])
  data.table(
    cell_id = meta$cell_id,
    state_final = meta$state_final,
    gene = gene_id,
    raw_count = raw
  )[, .(
    cells = .N,
    expressing_cells = sum(raw_count > 0),
    detection_fraction = mean(raw_count > 0),
    raw_transcripts = sum(raw_count),
    mean_raw_count = mean(raw_count)
  ), by = .(gene, state_final)]
}))
marker_metrics <- merge(marker_metrics, marker_table[, .(gene, marker_axis)], by = "gene")
fwrite(marker_metrics, file.path(table_dir, "dataset4_non_signature_marker_validation_detection_by_state.tsv"), sep = "\t")

theme_marker <- function() {
  theme_void(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 6)
    )
}

plot_marker <- function(gene_id, axis_label) {
  high_color <- if (grepl("^Intestinal", axis_label)) "#0a9396" else "#ae2012"
  FeaturePlot(
    obj,
    features = gene_id,
    reduction = "tsne",
    pt.size = 0.18,
    order = TRUE,
    min.cutoff = "q05",
    max.cutoff = "q95",
    cols = c("#eeeeee", high_color)
  ) +
    ggtitle(gene_id) +
    theme_marker()
}

plots <- lapply(seq_len(nrow(marker_table)), function(i) {
  plot_marker(marker_table$gene[i], marker_table$marker_axis[i])
})

pdf_path <- file.path(figure_dir, "GSE183904_non_signature_marker_validation_tsne.pdf")
pdf(pdf_path, width = 11, height = 5.6, onefile = TRUE, useDingbats = FALSE)
print(
  wrap_plots(plots, nrow = 2) +
    plot_annotation(
      title = "Dataset4 non-signature marker validation on final malignant epithelial t-SNE",
      subtitle = "Top row: intestinal markers not used for state assignment | Bottom row: diffuse/EMT markers not used for state assignment"
    )
)
dev.off()

message("Non-signature marker validation PDF written to: ", pdf_path)
message("Marker metrics written to: ", file.path(table_dir, "dataset4_non_signature_marker_validation_detection_by_state.tsv"))
