# Simple script to open the final GSE183904 object and plot the t-SNE.
# This only reads the RDS object. It does not modify the saved analysis.
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

rds_path <- "C:/Users/pcyou/Desktop/stage_LLM/M1_Bioinformatics_Ion_Channel/Gastric_KCNQ1/gastric_dataset4/outputs/objects/final_intestinal_diffuse_emt_seurat.rds"
obj <- readRDS(rds_path)

state_colors <- c(
  "Intestinal-like" = "#2b9348",
  "Diffuse/EMT-like" = "#ff5714"
)

# --- Compaction du t-SNE (uniquement pour l'affichage) ---
factor <- 0.5  # <1 = plus compact ; ajuste entre 0.3 (très compact) et 0.8 (léger)

coords <- Embeddings(obj, "tsne")
df <- data.frame(
  tSNE_1 = coords[, 1] * factor,
  tSNE_2 = coords[, 2] * factor,
  state  = obj$state_final
)

p_state <- ggplot(df, aes(x = tSNE_1, y = tSNE_2, color = state)) +
  geom_point(size = 0.4) +
  scale_color_manual(values = state_colors) +
  ggtitle("Data4 malignant epithelial states") +
  theme_void(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank()
  ) +
  coord_fixed()

print(p_state)




# Projection de l'expression de KCNQ1 sur le t-SNE
p_kcnq1 <- FeaturePlot(
  obj,
  features = "KCNQ1",
  reduction = "tsne",
  pt.size = 0.4
) +
  ggtitle("Data4 - KCNQ1 expression") +
  theme_void(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p_kcnq1)

# Optional: save
# ggsave("dataset4_KCNQ1_tsne.pdf", p_kcnq1, width = 7, height = 5)
