# Plot KCN expression gradients on the final dataset4 t-SNE.
# Yellow indicates low or no expression; purple indicates higher expression.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

root <- "C:/Users/pcyou/Desktop/stage_LLM/M1_Bioinformatics_Ion_Channel/Gastric_KCNQ1/gastric_dataset4"
rds_path <- file.path(root, "outputs/objects/final_intestinal_diffuse_emt_seurat.rds")
figure_dir <- file.path(root, "outputs/final/figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(rds_path)
DefaultAssay(obj) <- "RNA"

kcn_genes <- c("KCNQ1", "KCNE2", "KCNE3")
missing_genes <- setdiff(kcn_genes, rownames(obj))
if (length(missing_genes) > 0) {
  stop("Missing KCN gene(s): ", paste(missing_genes, collapse = ", "))
}

theme_kcn <- function() {
  theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 7)
    )
}

plots <- lapply(kcn_genes, function(gene_id) {
  FeaturePlot(
    obj,
    features = gene_id,
    reduction = "tsne",
    pt.size = 0.25,
    order = TRUE,
    min.cutoff = "q00",
    max.cutoff = "q95",
    cols = c("#ffd500", "#6247aa")
  ) +
    ggtitle(gene_id) +
    theme_kcn()
})

pdf_path <- file.path(figure_dir, "GSE183904_KCN_gradient_featureplots_tsne.pdf")
pdf(pdf_path, width = 9, height = 3.4, onefile = TRUE, useDingbats = FALSE)
print(
  wrap_plots(plots, nrow = 1) +
    plot_annotation(
      title = "Dataset4 KCN expression gradients on final malignant epithelial t-SNE",
      subtitle = "Yellow: low or no expression | Purple: higher log-normalized expression"
    )
)
dev.off()

message("KCN gradient FeaturePlot PDF written to: ", pdf_path)
