suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(infercnv)
})

set.seed(12345)
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
table_dir <- file.path(root, "outputs", "tables")
object_dir <- file.path(root, "outputs", "objects")
infer_dir <- file.path(root, "outputs", "intermediate", "infercnv")
dir.create(infer_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(file.path(object_dir, "epithelial_candidates_seurat.rds"))
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
gene_order_path <- file.path(table_dir, "infercnv_GRCh38_gene_order.tsv")
if (!file.exists(gene_order_path)) {
  stop("Missing GRCh38 gene-order file: ", gene_order_path,
       ". Keep this table in outputs/tables before rerunning inferCNV.")
}
gene_order <- fread(gene_order_path, header = FALSE, col.names = c("gene", "chr", "start", "end"))
gene_order <- unique(gene_order[gene %chin% rownames(counts)], by = "gene")
gene_order[, chr := sub("^chr", "", chr)]
gene_order <- gene_order[chr %chin% c(as.character(1:22), "X", "Y")]
gene_order[, chr_order := match(chr, c(as.character(1:22), "X", "Y"))]
setorder(gene_order, chr_order, start)
gene_order[, chr_order := NULL]

counts <- counts[gene_order$gene, , drop = FALSE]
keep_gene <- Matrix::rowSums(counts > 0) >= 20
counts <- counts[keep_gene, , drop = FALSE]
gene_order <- gene_order[gene %chin% rownames(counts)]
counts <- counts[gene_order$gene, , drop = FALSE]

meta <- as.data.table(obj@meta.data, keep.rownames = "cell_id")
reference_cells <- meta[tissue == "Primary_Normal", .SD[sample(.N, min(.N, 150L))], by = patient]$cell_id
observation_cells <- meta[tissue == "Primary_Tumor", cell_id]
analysis_cells <- c(reference_cells, observation_cells)
counts_infer <- counts[, analysis_cells, drop = FALSE]

annotations <- data.table(
  cell_id = analysis_cells,
  infercnv_group = c(rep("Reference_normal", length(reference_cells)), rep("Primary_tumor", length(observation_cells)))
)
annotation_path <- file.path(table_dir, "infercnv_annotations.tsv")
order_path <- file.path(table_dir, "infercnv_GRCh38_gene_order.tsv")
fwrite(annotations, annotation_path, sep = "\t", col.names = FALSE)
fwrite(gene_order, order_path, sep = "\t", col.names = FALSE)
fwrite(meta[cell_id %chin% reference_cells, .N, by = .(patient, lauren, tissue)],
       file.path(table_dir, "infercnv_reference_summary.tsv"), sep = "\t")

result_path <- file.path(object_dir, "infercnv_GSE183904.rds")
if (!file.exists(result_path)) {
  infer_obj <- CreateInfercnvObject(
    raw_counts_matrix = counts_infer,
    annotations_file = annotation_path,
    delim = "\t",
    gene_order_file = order_path,
    ref_group_names = "Reference_normal"
  )
  infer_res <- infercnv::run(
    infer_obj,
    cutoff = 0.1,
    out_dir = infer_dir,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = FALSE,
    num_threads = 8,
    no_plot = TRUE
  )
  saveRDS(infer_res, result_path, compress = FALSE)
} else {
  infer_res <- readRDS(result_path)
}

cnv_expr <- infer_res@expr.data
reference_present <- intersect(reference_cells, colnames(cnv_expr))
reference_center <- Matrix::rowMeans(cnv_expr[, reference_present, drop = FALSE])
cnv_deviation <- sweep(cnv_expr, 1, reference_center, FUN = "-")
cnv_score <- sqrt(colMeans(cnv_deviation^2))
cnv_table <- data.table(cell_id = names(cnv_score), cnv_score = as.numeric(cnv_score))
cnv_table <- merge(cnv_table, meta[, .(cell_id, sample_id, patient, tissue, lauren)],
                   by = "cell_id", all.x = TRUE, sort = FALSE)
reference_threshold <- quantile(cnv_table[tissue == "Primary_Normal", cnv_score], 0.95, na.rm = TRUE)
cnv_table[, cnv_high := cnv_score > reference_threshold]
fwrite(cnv_table, file.path(table_dir, "infercnv_cell_scores.tsv"), sep = "\t")
fwrite(data.table(
  parameter = c("reference_cells", "tumor_observation_cells", "genes_used", "normal_reference_q95"),
  value = c(length(reference_cells), length(observation_cells), nrow(counts_infer), reference_threshold)
), file.path(table_dir, "infercnv_run_summary.tsv"), sep = "\t")
message("inferCNV completed: ", length(reference_cells), " reference cells, ",
        length(observation_cells), " tumor cells and ", nrow(counts_infer), " genes.")
