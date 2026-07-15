library(Seurat)
library(data.table)
library(ggplot2)
library(patchwork)

base_dir <- "C:/Users/pcyou/Desktop/stage_LLM/M1_Bioinformatics_Ion_Channel/Gastric_KCNQ1/gastric_dataset4"
object_path <- file.path(base_dir, "outputs/objects/final_intestinal_diffuse_emt_seurat.rds")
out_dir <- file.path(base_dir, "outputs/final/figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(object_path)
emb <- as.data.table(Embeddings(obj, reduction = "tsne"), keep.rownames = "cell_id")
setnames(emb, old = c("tSNE_1", "tSNE_2"), new = c("tsne_1", "tsne_2"), skip_absent = TRUE)
if (!all(c("tsne_1", "tsne_2") %in% names(emb))) {
  coord_cols <- setdiff(names(emb), "cell_id")[1:2]
  setnames(emb, coord_cols, c("tsne_1", "tsne_2"))
}

meta <- as.data.table(obj@meta.data, keep.rownames = "cell_id")
dt <- merge(
  emb,
  meta[, .(cell_id, patient, lauren, state_final, seurat_clusters)],
  by = "cell_id",
  all.x = TRUE,
  sort = FALSE
)

dt[, state_final := factor(state_final, levels = c("Intestinal-like", "Diffuse/EMT-like"))]
dt[, patient := factor(patient)]
dt[, seurat_clusters := factor(seurat_clusters)]

cluster_centers <- dt[, .(
  tsne_1 = median(tsne_1),
  tsne_2 = median(tsne_2),
  n = .N,
  dominant_state = names(sort(table(state_final), decreasing = TRUE))[1],
  top_patient = names(sort(table(patient), decreasing = TRUE))[1]
), by = seurat_clusters]

theme_tsne_clean <- function() {
  theme_classic(base_size = 10) +
    theme(
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 7, color = "black"),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 8.5, hjust = 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = 5.5),
      legend.key.size = unit(0.24, "cm"),
      legend.spacing.y = unit(0.02, "cm"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

add_cluster_labels <- function(p) {
  p +
    geom_label(
      data = cluster_centers,
      aes(x = tsne_1, y = tsne_2, label = seurat_clusters),
      inherit.aes = FALSE,
      size = 2.6,
      label.size = 0.18,
      label.padding = unit(0.10, "lines"),
      fill = "white",
      alpha = 0.86,
      color = "black"
    )
}

state_cols <- c("Intestinal-like" = "#2b9348", "Diffuse/EMT-like" = "#ff5714")

p_state <- ggplot(dt, aes(tsne_1, tsne_2, color = state_final)) +
  geom_point(size = 0.35, alpha = 0.88, stroke = 0) +
  scale_color_manual(values = state_cols, drop = FALSE) +
  labs(
    title = "Final malignant epithelial states",
    subtitle = "High-confidence intestinal-like vs diffuse/EMT-like cells",
    x = "tSNE_1",
    y = "tSNE_2"
  ) +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1), ncol = 1)) +
  theme_tsne_clean() +
  coord_fixed(ratio = 0.72)
p_state <- add_cluster_labels(p_state)

p_patient <- ggplot(dt, aes(tsne_1, tsne_2, color = patient)) +
  geom_point(size = 0.30, alpha = 0.90, stroke = 0) +
  labs(
    title = "Same t-SNE colored by patient",
    subtitle = "Small legend to check patient-driven clusters",
    x = "tSNE_1",
    y = "tSNE_2"
  ) +
  guides(color = guide_legend(override.aes = list(size = 1.8, alpha = 1), ncol = 2)) +
  theme_tsne_clean() +
  coord_fixed(ratio = 0.72) +
  theme(legend.position = "right")
p_patient <- add_cluster_labels(p_patient)

p_cluster <- ggplot(dt, aes(tsne_1, tsne_2, color = seurat_clusters)) +
  geom_point(size = 0.32, alpha = 0.86, stroke = 0) +
  labs(
    title = "Unsupervised clusters",
    subtitle = "Cluster numbers shown at median t-SNE position",
    x = "tSNE_1",
    y = "tSNE_2"
  ) +
  guides(color = "none") +
  theme_tsne_clean() +
  coord_fixed(ratio = 0.72)
p_cluster <- add_cluster_labels(p_cluster)

pdf_path <- file.path(out_dir, "GSE183904_final_states_patient_cluster_labels_tsne.pdf")
pdf(pdf_path, width = 12.6, height = 6.2, useDingbats = FALSE)
print(
  (p_state | p_patient | p_cluster) +
    plot_annotation(
      title = "Dataset4 final malignant-cell states: state, patient and cluster structure",
      theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
    )
)
dev.off()

fwrite(
  cluster_centers[order(as.integer(as.character(seurat_clusters)))],
  file.path(out_dir, "GSE183904_final_state_cluster_label_summary.tsv"),
  sep = "\t"
)

message("Saved: ", pdf_path)
