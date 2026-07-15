suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

set.seed(12345)
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
table_dir <- file.path(root, "outputs", "tables")
object_dir <- file.path(root, "outputs", "objects")
figure_dir <- file.path(root, "outputs", "intermediate", "figures")
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

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

object_path <- file.path(object_dir, "epithelial_candidates_seurat.rds")
if (file.exists(object_path)) {
  obj <- readRDS(object_path)
} else {
  counts <- readMM(file.path(root, "outputs", "intermediate", "epithelial_counts.mtx"))
  features <- fread(file.path(table_dir, "features.tsv"))$gene
  manifest <- fread(file.path(table_dir, "epithelial_cell_manifest.tsv"))
  stopifnot(nrow(counts) == length(features), ncol(counts) == nrow(manifest))
  rownames(counts) <- make.unique(features)
  colnames(counts) <- manifest$cell_id
  counts <- as(counts, "dgCMatrix")

  obj <- CreateSeuratObject(counts = counts, project = "GSE183904_epithelial", min.cells = 3)
  meta <- as.data.frame(manifest)
  rownames(meta) <- meta$cell_id
  obj <- AddMetaData(obj, meta[colnames(obj), setdiff(names(meta), "cell_id"), drop = FALSE])
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, vars.to.regress = "nCount_RNA", features = VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30, seed.use = 12345, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:20, k.param = 30, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.8, random.seed = 12345, verbose = FALSE)
  obj <- RunTSNE(obj, dims = 1:20, seed.use = 12345, check_duplicates = FALSE, verbose = FALSE)
  saveRDS(obj, object_path, compress = FALSE)
}

cluster_summary <- as.data.table(obj@meta.data, keep.rownames = "cell_id")[, .(
  cells = .N,
  tumor_fraction = mean(tissue == "Primary_Tumor"),
  dominant_patient_fraction = max(table(patient)) / .N,
  median_features = as.numeric(median(nFeature_RNA)),
  median_counts = as.numeric(median(nCount_RNA))
), by = seurat_clusters]
fwrite(cluster_summary, file.path(table_dir, "epithelial_cluster_summary.tsv"), sep = "\t")

p1 <- DimPlot(obj, reduction = "tsne", group.by = "tissue", pt.size = 0.16,
              cols = c("Primary_Normal" = "#88d498", "Primary_Tumor" = "#ff6b6b")) +
  labs(title = "Broad epithelial candidates", subtitle = "Primary tumor and matched normal") + theme_tsne()
p2 <- DimPlot(obj, reduction = "tsne", group.by = "lauren", pt.size = 0.16,
              cols = c("Intestinal" = "#0a9396", "Diffuse" = "#ae2012")) +
  labs(title = "Patient Lauren type", subtitle = "Clinical metadata, not the cell-state classifier") + theme_tsne()
p3 <- DimPlot(obj, reduction = "tsne", group.by = "patient", pt.size = 0.12) +
  labs(title = "Patient origin", subtitle = "Interpatient structure is retained") + theme_tsne()
p4 <- DimPlot(obj, reduction = "tsne", group.by = "seurat_clusters", pt.size = 0.16,
              label = TRUE, repel = TRUE) +
  labs(title = "Unsupervised epithelial clusters") + theme_tsne()
page <- wrap_plots(lapply(list(p1, p2, p3, p4), single_plot), ncol = 2) +
  plot_annotation(title = "GSE183904 epithelial-cell preparation")
pdf(file.path(figure_dir, "01_epithelial_preparation_tsne.pdf"), width = 10, height = 8,
    onefile = TRUE, useDingbats = FALSE)
print(page)
dev.off()
message("Epithelial object completed: ", ncol(obj), " cells and ", nrow(obj), " retained genes.")
