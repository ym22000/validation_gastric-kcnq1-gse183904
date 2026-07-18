suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
analysis_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
dataset_root <- normalizePath(file.path(analysis_root, ".."), winslash = "/")
figure_dir <- file.path(analysis_root, "outputs", "figures")
table_dir <- file.path(analysis_root, "outputs", "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(file.path(dataset_root, "outputs", "objects", "final_intestinal_diffuse_emt_seurat.rds"))
DefaultAssay(obj) <- "RNA"
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
expression <- GetAssayData(obj, assay = "RNA", layer = "data")
coordinates <- as.data.table(Embeddings(obj, "tsne"), keep.rownames = "cell_id")
setnames(coordinates, c("cell_id", "x", "y"))

genes <- c("KCNQ1", "KCNE3", "KCNQ2")
stopifnot(all(genes %in% rownames(counts)))

membership <- fread(file.path(table_dir, "cohort_cell_membership.tsv.gz"))
cohort_order <- c("Original final", "Cluster-screened", "Orthogonal strict")
membership[, cohort := factor(cohort, levels = cohort_order)]

compact_coordinates <- function(data) {
  data[, x := x - mean(range(x))]
  data[, y := y - mean(range(y))]
  common_scale <- max(diff(range(data$x)), diff(range(data$y)))
  data[, `:=`(x = x / common_scale, y = y / common_scale)]
  data
}

plot_gene_cohort <- function(gene, cohort_name, upper_limit) {
  ids <- membership[cohort == cohort_name, cell_id]
  data <- coordinates[cell_id %chin% ids]
  data[, expression := as.numeric(expression[gene, cell_id])]
  data[, raw_count := as.numeric(counts[gene, cell_id])]
  data <- compact_coordinates(data)
  setorder(data, expression)
  detection <- mean(data$raw_count > 0)

  ggplot(data, aes(x, y, color = expression)) +
    geom_point(size = 0.34, alpha = 0.92) +
    scale_color_gradient(
      low = "#ffd500", high = "#6247aa",
      limits = c(0, upper_limit), oob = scales::squish,
      name = "Log-normalized\nexpression"
    ) +
    coord_fixed(xlim = c(-0.54, 0.54), ylim = c(-0.54, 0.54), clip = "off") +
    labs(
      title = cohort_name,
      subtitle = sprintf("%s cells | %.2f%% detected", format(nrow(data), big.mark = ","), 100 * detection)
    ) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = element_text(size = 8.5, hjust = 0.5, color = "grey30"),
      plot.margin = margin(4, 5, 4, 5),
      legend.title = element_text(face = "bold", size = 8.5),
      legend.text = element_text(size = 8)
    )
}

pdf_path <- file.path(figure_dir, "GSE183904_KCNQ1_KCNE3_KCNQ2_three_cohort_tsne.pdf")
pdf(pdf_path, width = 9.6, height = 3.8, onefile = TRUE, useDingbats = FALSE)
for (gene in genes) {
  positive_values <- as.numeric(expression[gene, colnames(obj)])
  positive_values <- positive_values[positive_values > 0]
  upper_limit <- if (length(positive_values) >= 10L) {
    as.numeric(quantile(positive_values, 0.99, na.rm = TRUE))
  } else if (length(positive_values)) {
    max(positive_values)
  } else {
    1
  }
  upper_limit <- max(upper_limit, .Machine$double.eps)

  page <- wrap_plots(
    lapply(cohort_order, function(cohort) plot_gene_cohort(gene, cohort, upper_limit)),
    nrow = 1,
    guides = "collect"
  ) +
    plot_annotation(title = gene) &
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      legend.position = "right"
    )
  print(page)
}
dev.off()

projection_summary <- rbindlist(lapply(genes, function(gene) {
  membership[, .(
    cells = .N,
    expressing_cells = sum(as.numeric(counts[gene, cell_id]) > 0),
    detection_fraction = mean(as.numeric(counts[gene, cell_id]) > 0),
    raw_transcripts = sum(as.numeric(counts[gene, cell_id]))
  ), by = cohort][, gene := gene]
}))
setcolorder(projection_summary, c("gene", "cohort", "cells", "expressing_cells", "detection_fraction", "raw_transcripts"))
fwrite(projection_summary, file.path(table_dir, "KCN_three_cohort_projection_summary.tsv"), sep = "\t")

message("Three-cohort KCN t-SNE PDF written to: ", pdf_path)

