# Merge E1 + E2 + E3 annotated Seurat objects, SCTransform, and Harmony integration.
# Run with working directory set to the project root.
# Prerequisites: Barcode_E1_Clones.Rmd, Barcode_E2_Clones.Rmd, and Barcode_E3_Clones.Rmd
# (annotated objects with clone metadata in Data/E16_barcodes_E{1,2,3}.rds).

library(Seurat)
options(Seurat.object.assay.version = "v5")
library(sctransform)
library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)
library(reticulate)
library(cowplot)
library(scCustomize)
py_module_available(module = "leidenalg")
options(future.globals.maxSize = 100 * 1024^3)

fig_dir_regions <- "Figures/Fig5/"

save_fig_regions <- function(name, width = 7, height = 5, device = grDevices::pdf, bg = "white", dpi = 300) {
  if (!dir.exists(fig_dir_regions)) dir.create(fig_dir_regions, recursive = TRUE)

  ggsave(
    filename = file.path(fig_dir_regions, name),
    width = width,
    height = height,
    bg = bg,
    dpi = dpi,
    device = device
  )
}

for (f in c(
  "Data/E16_barcodes_E1.rds",
  "Data/E16_barcodes_E2.rds",
  "Data/E16_barcodes_E3.rds"
)) {
  if (!file.exists(f)) {
    stop("Missing ", f, ". Run the corresponding processing and clonal notebooks first.")
  }
}

# Importing all data ---------------------------------------------------------------
E1 <- readRDS("Data/E16_barcodes_E1.rds")
E2 <- readRDS("Data/E16_barcodes_E2.rds")
E3 <- readRDS("Data/E16_barcodes_E3.rds")

Idents(E1) <- "assigned_cell_type"
Idents(E2) <- "assigned_cell_type"
Idents(E3) <- "assigned_cell_type"

DefaultAssay(E1) <- "RNA"
DefaultAssay(E2) <- "RNA"
DefaultAssay(E3) <- "RNA"

E1[["SCT"]] <- NULL
E2[["SCT"]] <- NULL
E3[["SCT"]] <- NULL

Barcodes_combined <- merge(E1, y = c(E2, E3), add.cell.ids = c("E1", "E2", "E3"), project = "Barcoded")
Barcodes_combined$embryo <- sub("_.*", "", colnames(Barcodes_combined))

Barcodes_combined <- SCTransform(Barcodes_combined, vars.to.regress = "percent.mt", verbose = TRUE)
VarGenes_Barcodes_combined <- VariableFeatures(object = Barcodes_combined)
VarGenes_Barcodes_combined <- VarGenes_Barcodes_combined[!VarGenes_Barcodes_combined %in% c(
  "Xist", "Tsix", "Ddx3y", "Uty", "Eif2s3y", "Kdm5d", "Gm29650", "Malat1",
  "Gm42418", "AY036118", "Gm13305"
)]
Barcodes_combined <- RunPCA(Barcodes_combined, verbose = FALSE, features = VarGenes_Barcodes_combined)
Barcodes_combined <- RunUMAP(
  Barcodes_combined, reduction = "pca", dims = 1:50,
  n.neighbors = 55L, min.dist = 0.40, n.epochs = 1000, local.connectivity = 1L, seed.use = 11
)
DimPlot(Barcodes_combined, label = TRUE, group.by = "embryo", reduction = "umap")
DimPlot(Barcodes_combined, label = TRUE, group.by = "sample", reduction = "umap")

Barcodes_combined <- IntegrateLayers(
  object = Barcodes_combined, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  assay = "SCT", verbose = FALSE
)

Barcodes_combined <- RunUMAP(
  Barcodes_combined, reduction = "harmony", dims = 1:50, reduction.name = "umap.harmony",
  n.neighbors = 55L, min.dist = 0.40, n.epochs = 1000, local.connectivity = 1L
)
Barcodes_combined <- FindNeighbors(Barcodes_combined, reduction = "harmony", dims = 1:50)
Barcodes_combined <- FindClusters(
  Barcodes_combined, resolution = 1.7, cluster.name = "harmony_clusters",
  algorithm = 4, random.seed = 111
)
Barcodes_combined$embryo <- factor(
  Barcodes_combined$embryo,
  levels = c("E3", "E2", "E1")
)
DimPlot(Barcodes_combined, reduction = "umap.harmony", group.by = "embryo")
save_fig_regions("Merged_byregion_3datasets.pdf", width = 8, height = 6)

DimPlot(Barcodes_combined, reduction = "umap.harmony", group.by = "Phase")
DimPlot(Barcodes_combined, reduction = "umap.harmony", group.by = "clones")
DimPlot(Barcodes_combined, reduction = "umap.harmony", group.by = "sample")
DimPlot(Barcodes_combined, reduction = "umap.harmony", label = TRUE)
DimPlot(Barcodes_combined, reduction = "umap.harmony", label = TRUE, group.by = "major_cell_type")
DimPlot(Barcodes_combined, label = TRUE, group.by = "assigned_cell_type", reduction = "umap.harmony")
VlnPlot(Barcodes_combined, features = c("nUMI", "nGene", "percent.mt"), ncol = 3, pt.size = 0.01)

saveRDS(object = Barcodes_combined, "Data/E16_barcodes_E1E2E3_merged.rds")

Barcodes_combined <- PrepSCTFindMarkers(Barcodes_combined)
Barcodes_combined_markers <- FindAllMarkers(
  Barcodes_combined, only.pos = TRUE, test.use = "wilcox",
  min.pct = 0.25, logfc.threshold = 0.25, group.by = "harmony_clusters"
)

top30_Barcodes_combined_markers <- Barcodes_combined_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 30)

write.csv(
  top30_Barcodes_combined_markers,
  file = "Data/top30_Barcodes_combined_markers.csv",
  row.names = FALSE
)
