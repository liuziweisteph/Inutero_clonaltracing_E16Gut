library(Seurat)

# Input and output files
input_file <- "E16_barcodes_E2_202512.rds"
output_file <- "seurat_E16_barcodes_E2_202512_harmony_clusters.csv"

cat("Reading RDS file:", input_file, "\n")
seurat_obj <- readRDS(input_file)

cat("Loaded:", ncol(seurat_obj), "cells,", nrow(seurat_obj), "genes\n")

DefaultAssay(seurat_obj) <- DefaultAssay(seurat_obj)
if ("SCT" %in% names(seurat_obj@assays)) {
    cat("SCT assay detected. Running PrepSCTFindMarkers...\n")
    seurat_obj <- PrepSCTFindMarkers(seurat_obj)
}

# Set identity to cluster column (using group.by in FindAllMarkers)
cat("Running FindAllMarkers with group.by = 'harmony_clusters'...\n")


markers <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    test.use = "wilcox",
    min.pct = 0.25,
    logfc.threshold = 0.5,
    group.by = "harmony_clusters"
)


cat("Columns in markers dataframe:", paste(colnames(markers), collapse = ", "), "\n")
cat("Number of markers:", nrow(markers), "\n")

output_df <- data.frame(
    p_val = markers$p_val,
    avg_logFC = if("avg_log2FC" %in% colnames(markers)) markers$avg_log2FC else markers$avg_logFC,
    pct.1 = markers$pct.1,
    pct.2 = markers$pct.2,
    p_val_adj = markers$p_val_adj,
    cluster = markers$cluster,
    gene = markers$gene,
    stringsAsFactors = FALSE
)

# Sort by cluster, then p_val
output_df <- output_df[order(output_df$cluster, output_df$p_val), ]
write.csv(output_df, file = output_file, row.names = FALSE, quote = FALSE)
cat("Done! Saved", nrow(output_df), "markers to", output_file, "\n")


