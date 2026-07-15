library(Seurat)
library(data.table)
library(msigdbr)
library(fgsea)
library(ggplot2)

base_dir <- "C:/Users/pcyou/Desktop/stage_LLM/M1_Bioinformatics_Ion_Channel/Gastric_KCNQ1/gastric_dataset4"
object_path <- file.path(base_dir, "outputs/objects/final_intestinal_diffuse_emt_seurat.rds")
out_dir <- file.path(base_dir, "outputs/final/gsea_hallmark_selected_clusters")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

diffuse_clusters <- c("2", "5", "17", "18", "21")
intestinal_clusters <- c("7", "11", "12", "8", "0", "15", "14")

obj <- readRDS(object_path)
meta <- as.data.table(obj@meta.data, keep.rownames = "cell_id")
meta[, seurat_clusters := as.character(seurat_clusters)]
meta[, selected_pole := fifelse(
  seurat_clusters %in% diffuse_clusters,
  "Diffuse_EMT_selected_clusters",
  fifelse(seurat_clusters %in% intestinal_clusters, "Intestinal_selected_clusters", "Other_final_clusters")
)]

cluster_summary <- meta[, .(
  cells = .N,
  intestinal_like = sum(state_final == "Intestinal-like"),
  diffuse_emt_like = sum(state_final == "Diffuse/EMT-like"),
  n_patients = uniqueN(patient),
  top_patient = names(sort(table(patient), decreasing = TRUE))[1],
  top_patient_fraction = as.numeric(max(table(patient)) / .N),
  lauren_breakdown = paste(names(table(lauren)), as.integer(table(lauren)), sep = ":", collapse = ";")
), by = .(seurat_clusters, selected_pole)][order(selected_pole, as.integer(seurat_clusters))]

fwrite(cluster_summary, file.path(out_dir, "selected_cluster_composition.tsv"), sep = "\t")

selected_cells <- meta[selected_pole != "Other_final_clusters", cell_id]
obj_sel <- subset(obj, cells = selected_cells)
obj_sel$selected_pole <- meta[match(colnames(obj_sel), cell_id), selected_pole]
Idents(obj_sel) <- "selected_pole"

deg <- FindMarkers(
  object = obj_sel,
  ident.1 = "Diffuse_EMT_selected_clusters",
  ident.2 = "Intestinal_selected_clusters",
  assay = "RNA",
  slot = "data",
  test.use = "wilcox",
  min.pct = 0,
  logfc.threshold = 0,
  only.pos = FALSE,
  verbose = FALSE
)
deg <- as.data.table(deg, keep.rownames = "gene")

fc_col <- intersect(c("avg_log2FC", "avg_logFC"), names(deg))[1]
if (is.na(fc_col)) stop("No average logFC column found in FindMarkers output.")
setnames(deg, fc_col, "avg_log2FC")
deg[, comparison := "Diffuse_EMT_selected_clusters_vs_Intestinal_selected_clusters"]
deg[, direction := fifelse(avg_log2FC > 0, "Higher in selected Diffuse/EMT clusters",
                           fifelse(avg_log2FC < 0, "Higher in selected Intestinal clusters", "No direction"))]
deg[, rank_metric := sign(avg_log2FC) * -log10(pmax(p_val, 1e-300))]
setorder(deg, -rank_metric)
fwrite(deg, file.path(out_dir, "DEG_selected_diffuse_EMT_vs_intestinal_clusters.tsv"), sep = "\t")

rank_dt <- deg[is.finite(rank_metric) & !is.na(gene) & gene != ""]
rank_dt <- rank_dt[order(abs(rank_metric), decreasing = TRUE)]
rank_dt <- rank_dt[!duplicated(gene)]
ranks <- rank_dt$rank_metric
names(ranks) <- rank_dt$gene
ranks <- sort(ranks, decreasing = TRUE)

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)
gene_sets_used <- data.table(
  pathway = names(pathways),
  n_genes_total = lengths(pathways),
  n_genes_in_rank = vapply(pathways, function(g) sum(g %in% names(ranks)), integer(1))
)
fwrite(gene_sets_used, file.path(out_dir, "hallmark_gene_sets_used.tsv"), sep = "\t")

fg <- fgsea(
  pathways = pathways,
  stats = ranks,
  minSize = 15,
  maxSize = 500,
  eps = 0
)
fg <- as.data.table(fg)
fg[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
fg[, enriched_in := fifelse(NES > 0, "Selected Diffuse/EMT clusters", "Selected Intestinal clusters")]
fg[, abs_NES := abs(NES)]
setorder(fg, padj, -abs_NES)
fwrite(fg, file.path(out_dir, "hallmark_fgsea_results.tsv"), sep = "\t")

top_plot <- rbind(
  fg[NES > 0][order(padj, -NES)][1:min(.N, 12)],
  fg[NES < 0][order(padj, NES)][1:min(.N, 12)]
)
top_plot <- top_plot[!is.na(pathway)]
top_plot[, pathway_clean := gsub("^HALLMARK_", "", pathway)]
top_plot[, pathway_clean := gsub("_", " ", pathway_clean)]
top_plot[, pathway_clean := fifelse(padj < 0.05, paste0(pathway_clean, " *"), pathway_clean)]
top_plot[, pathway_clean := factor(pathway_clean, levels = pathway_clean[order(NES)])]

pdf(file.path(out_dir, "hallmark_gsea_selected_clusters_summary.pdf"), width = 8.4, height = 7.2, useDingbats = FALSE)
print(
  ggplot(top_plot, aes(x = pathway_clean, y = NES, fill = enriched_in)) +
    geom_col(width = 0.75) +
    coord_flip() +
    scale_fill_manual(values = c(
      "Selected Diffuse/EMT clusters" = "#AE2012",
      "Selected Intestinal clusters" = "#0A9396"
    )) +
    labs(
      title = "Hallmark GSEA on selected final malignant-cell clusters",
      subtitle = "Positive NES: selected diffuse/EMT clusters | Negative NES: selected intestinal clusters",
      x = NULL,
      y = "Normalized enrichment score (NES)",
      fill = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom",
      axis.text.y = element_text(size = 7.5)
    )
)
dev.off()

top_deg_diffuse <- deg[avg_log2FC > 0][order(p_val_adj, -avg_log2FC)][1:50]
top_deg_intestinal <- deg[avg_log2FC < 0][order(p_val_adj, avg_log2FC)][1:50]
fwrite(top_deg_diffuse, file.path(out_dir, "top50_genes_higher_in_selected_diffuse_EMT_clusters.tsv"), sep = "\t")
fwrite(top_deg_intestinal, file.path(out_dir, "top50_genes_higher_in_selected_intestinal_clusters.tsv"), sep = "\t")

summary_txt <- file.path(out_dir, "analysis_summary.txt")
writeLines(c(
  "Hallmark GSEA on selected dataset4 final malignant-cell clusters",
  "",
  paste0("Diffuse/EMT selected clusters: ", paste(diffuse_clusters, collapse = ", ")),
  paste0("Intestinal selected clusters: ", paste(intestinal_clusters, collapse = ", ")),
  paste0("Cells in selected diffuse/EMT clusters: ", ncol(obj_sel[, obj_sel$selected_pole == "Diffuse_EMT_selected_clusters"])),
  paste0("Cells in selected intestinal clusters: ", ncol(obj_sel[, obj_sel$selected_pole == "Intestinal_selected_clusters"])),
  "",
  "Differential expression:",
  "Seurat FindMarkers, Wilcoxon rank-sum test, RNA log-normalized data.",
  "Positive avg_log2FC means higher in selected diffuse/EMT clusters.",
  "Negative avg_log2FC means higher in selected intestinal clusters.",
  "",
  "GSEA:",
  "Hallmark gene sets from msigdbr, Homo sapiens, collection H.",
  "fgsea minSize=15, maxSize=500, eps=0.",
  "Ranking metric = sign(avg_log2FC) * -log10(raw p-value capped at 1e-300).",
  "Positive NES means enrichment in selected diffuse/EMT clusters.",
  "Negative NES means enrichment in selected intestinal clusters."
), summary_txt)

message("Saved outputs in: ", out_dir)
