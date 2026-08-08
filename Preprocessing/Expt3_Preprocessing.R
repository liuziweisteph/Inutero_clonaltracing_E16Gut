# Preprocessing Experiment 3: SoupX correction, QC, doublet removal, merge, SCTransform, Harmony, label transfer
# Run this script with working directory set to the project root

library(Seurat)
options(Seurat.object.assay.version = "v5")
library(DoubletFinder)
library(sctransform)
library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)
library(reticulate)
library(cowplot)
library(SoupX)
library(DropletUtils)
library(scCustomize)
py_module_available(module = "leidenalg")
options(future.globals.maxSize = 100 * 1024^3)
source("Helpers/SoupXDoubletfinder.R")

# External data directory — edit this path to where you extracted the raw Cell Ranger outputs (from GEO GSE325733)
external_data_dir <- "~/scRNA_seq/Trex_ENS/Submission/Data"
fig_dir_E3 <- "Figures/Fig6/"

save_fig_E3 <- function(name, width = 7, height = 5, device = grDevices::pdf, bg = "white", dpi = 300) {
  if (!dir.exists(fig_dir_E3)) dir.create(fig_dir_E3, recursive = TRUE)
  ggsave(
    filename = file.path(fig_dir_E3, name),
    width = width,
    height = height,
    bg = bg,
    dpi = dpi,
    device = device
  )
}

# SoupX -------------------------------------------------------------------
in_data_dir <- file.path(external_data_dir, "counts/E3")
soupx_out_dir <- file.path(external_data_dir, "SoupX/E3")

soupX_correction_counts_andsave(in_data_dir, soupx_out_dir)

# Experiment 3 ---------------------------------------------------------------
X10_e16_s1_E3.data <- Read10X(data.dir = file.path(soupx_out_dir, "E16_stomach_E3_SoupX_corrected"))
e16_s1_E3 <- CreateSeuratObject(counts = X10_e16_s1_E3.data, min.cells = 3, min.features = 200, project = "E16_stomach_barcodes_E3")
e16_s1_E3 <- PercentageFeatureSet(e16_s1_E3, pattern = "^mt-", col.name = "percent.mt")
VlnPlot(object = e16_s1_E3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
e16_s1_E3 <- subset(e16_s1_E3, subset = nFeature_RNA > 2500 & nFeature_RNA < 5000 & percent.mt < 10 & nCount_RNA > 2000 & nCount_RNA < 15000) #6174 cells
VlnPlot(object = e16_s1_E3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
rm(X10_e16_s1_E3.data)

X10_e16_s2_E3.data <- Read10X(data.dir = file.path(soupx_out_dir, "E16_sm_E3_SoupX_corrected"))
e16_s2_E3 <- CreateSeuratObject(counts = X10_e16_s2_E3.data, min.cells = 3, min.features = 200, project = "E16_jejileum_barcodes_E3")
e16_s2_E3 <- PercentageFeatureSet(e16_s2_E3, pattern = "^mt-", col.name = "percent.mt")
VlnPlot(object = e16_s2_E3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
e16_s2_E3 <- subset(e16_s2_E3, subset = nFeature_RNA > 2500 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA > 2000 & nCount_RNA < 17000) # 9563 cells
VlnPlot(object = e16_s2_E3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
rm(X10_e16_s2_E3.data)

X10_e16_s3_E3.data <- Read10X(data.dir = file.path(soupx_out_dir, "E16_colon_E3_SoupX_corrected"))
e16_s3_E3 <- CreateSeuratObject(counts = X10_e16_s3_E3.data, min.cells = 3, min.features = 200, project = "E16_colon_barcodes_E3")
e16_s3_E3 <- PercentageFeatureSet(e16_s3_E3, pattern = "^mt-", col.name = "percent.mt")
VlnPlot(object = e16_s3_E3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
e16_s3_E3 <- subset(e16_s3_E3, subset = nFeature_RNA > 2500 & nFeature_RNA < 7500 & percent.mt < 10 & nCount_RNA > 2000 & nCount_RNA < 50000) # 3260 cells
VlnPlot(object = e16_s3_E3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
rm(X10_e16_s3_E3.data)

E16_E3_list <- list("Stomach_E3" = e16_s1_E3, "Jej_ile_E3" = e16_s2_E3, "colon_E3" = e16_s3_E3)

# Doublet detection per sample ------------------------------------------------------
expectRate_E3 <- c(0.048, 0.072, 0.024)
names(expectRate_E3) <- c("Stomach_E3", "Jej_ile_E3", "colon_E3")

E16_E3_list <- process_seurat_objects_fordb(E16_E3_list)
E16_E3_list <- do_doubletfinder(E16_E3_list, expectRate_E3)

# Remove doublets ---------------------------------------------------------
E16_E3_list_clean <- sapply(E16_E3_list, function(seurat_obj) {
  subset(seurat_obj, subset = db_class == "Singlet")
}, simplify = FALSE)

cell_counts_E3 <- sapply(E16_E3_list_clean, function(seurat_obj) {
  ncol(seurat_obj)
})

print(cell_counts_E3) # After SoupX Stomach_E3 5904 Jej_ile_E3 8930 colon_E3 3189
# Merge datasets  ---------------------------------------------------------
E16_barcodes_E3 <- merge(E16_E3_list_clean[[1]], y = c(E16_E3_list_clean[[2]], E16_E3_list_clean[[3]]), add.cell.ids = c("e16_s1", "e16_ji1", "e16_c1"), project = "ENS_E16_barcodes_E3")

E16_barcodes_metadata_E3 <- E16_barcodes_E3@meta.data
# Add cell IDs to metadata
E16_barcodes_metadata_E3$cells <- rownames(E16_barcodes_metadata_E3)

# Rename columns
E16_barcodes_metadata_E3 <- E16_barcodes_metadata_E3 %>%
  dplyr::rename(
    seq_folder = orig.ident,
    nUMI = nCount_RNA,
    nGene = nFeature_RNA
  )

# Create sample column
E16_barcodes_metadata_E3$sample <- NA
E16_barcodes_metadata_E3$sample[which(str_detect(E16_barcodes_metadata_E3$cells, "^e16_s1_"))] <- "E16_stomach"
E16_barcodes_metadata_E3$sample[which(str_detect(E16_barcodes_metadata_E3$cells, "^e16_ji1_"))] <- "E16_jej_ileum"
E16_barcodes_metadata_E3$sample[which(str_detect(E16_barcodes_metadata_E3$cells, "^e16_c1_"))] <- "E16_colon"
E16_barcodes_metadata_E3$mitoRatio <- E16_barcodes_metadata_E3$percent.mt / 100
E16_barcodes_metadata_E3$log10GenesPerUMI <- log10(E16_barcodes_metadata_E3$nGene) / log10(E16_barcodes_metadata_E3$nUMI)
E16_barcodes_metadata_E3$batch <- "Embryo 3"
E16_barcodes_E3@meta.data <- E16_barcodes_metadata_E3

E16_barcodes_E3 <- SCTransform(E16_barcodes_E3, vars.to.regress = c("percent.mt", "nGene","nUMI"), verbose = TRUE)
VarGenes_E16_barcodes_E3 <- VariableFeatures(object = E16_barcodes_E3)
VarGenes_E16_barcodes_E3 <- VarGenes_E16_barcodes_E3[!VarGenes_E16_barcodes_E3 %in% c("Xist", "Gm13305", "Tsix", "Gm8730", "Eif2s3y", "Ddx3y", "Uty", "Kdm5d")]
E16_barcodes_E3 <- RunPCA(E16_barcodes_E3, verbose = FALSE, features = VarGenes_E16_barcodes_E3)
E16_barcodes_E3 <- RunUMAP(E16_barcodes_E3, reduction = "pca", dims = 1:50, n.neighbors = 55L, min.dist = 0.40, n.epochs = 1000, local.connectivity = 1L, seed.use = 111)
DimPlot(E16_barcodes_E3, reduction = "umap", group.by = "sample")

# Integration using Harmony
E16_barcodes_E3 <- IntegrateLayers(
  object = E16_barcodes_E3, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  assay = "SCT", verbose = FALSE, group.by.vars = "sample"
)
ElbowPlot_harmony(E16_barcodes_E3, ndims = 50, reduction = "harmony") # theta = 3
E16_barcodes_E3 <- RunUMAP(E16_barcodes_E3, reduction = "harmony", dims = 1:50, reduction.name = "umap.harmony", random.seed = 222,n.neighbors = 40L, min.dist = 0.10, n.epochs = 1000) # set neighbour =30, SCP disppeared
E16_barcodes_E3 <- FindNeighbors(E16_barcodes_E3, reduction = "harmony", dims = 1:50) # set neighbour =30, SCP disappeared
E16_barcodes_E3 <- FindClusters(E16_barcodes_E3, resolution = 1.7, cluster.name = "harmony_clusters", algorithm = 4, random.seed = 222)
DimPlot(E16_barcodes_E3, reduction = "umap.harmony", group.by = "sample")
DimPlot(E16_barcodes_E3, reduction = "umap.harmony", label = TRUE)

VlnPlot(E16_barcodes_E3,features = c("nUMI", "nGene", "percent.mt"), ncol = 3, pt.size = 0.01)

# Cell cycle scoring (gene lists from Expt1_Preprocessing.R → Data/mmus_*_genes.txt)
for (f in c("Data/mmus_s_genes.txt", "Data/mmus_g2m_genes.txt")) {
  if (!file.exists(f)) {
    stop("Missing ", f, ". Run Expt1_Preprocessing.R first to generate cell-cycle gene lists.")
  }
}
mmus_s <- read.table("Data/mmus_s_genes.txt",
  sep = "\t", header = FALSE, stringsAsFactors = FALSE
)$V1
mmus_g2m <- read.table("Data/mmus_g2m_genes.txt",
  sep = "\t", header = FALSE, stringsAsFactors = FALSE
)$V1
mmus_s <- mmus_s[-1]
mmus_g2m <- mmus_g2m[-1]

E16_barcodes_E3 <- CellCycleScoring(E16_barcodes_E3,
  s.features = mmus_s,
  g2m.features = mmus_g2m, set.ident = TRUE
)

DimPlot(E16_barcodes_E3, group.by = "Phase", reduction = "umap.harmony")
Idents(E16_barcodes_E3) <- E16_barcodes_E3$harmony_clusters

cell_type_features <- list(
  Muscle = c("Myh11", "Acta2", "Actg2"),
  Fibroblasts = c("Col1a1", "Lum"),
  Epithelial = c("Cdh1", "Epcam"),
  Endothelial = c("Cd36", "Pecam1"),
  Mesothelial = c("Upk3b", "Msln", "Wt1"),
  ENS = c("Phox2b", "Hand2", "Ret"),
  ICC = c("Kit", "Ano1"),
  Immnue = c("Lyz2", "Mrc1", "Trbc2", "Cd52", "Ptprc"),
  Pericytes = c("Rgs5", "Reln", "Pdgfrb", "Pdgfra"), # Pdgfra negative
  Myofibroblast = c("Trpc3", "Acta2", "Tagln"),
  Proliferation = c("Top2a", "Mki67"),
  Enterocytes = c("Vil1", "Krt20"),
  Enteroendocrine = c("Chga", "Tph1"),
  Goblet = c("Fcgbp", "Clca1"),
  Stomach = c("Pitx1")
)

ENS_features <- list(
  ENS = c("Phox2b", "Hand2", "Ret"),   
  Glia = c("Gfap", "Hey2", "Entpd2", "Heyl", "Eln", "Apoe"),
  Progenitor = c("Foxd3", "Erbb3", "Sox10", "Top2a", "Olfml3", "Col18a1", "Metrn"),
  SCP = c("Dhh", "Col14a1", "Gfra3", "Mal"),
  Neuroblast = c("Ascl1", "Insm1", "Dll3", "Bcl11b", "Hes6", "Btbd17"),
  Neuron = c("Stmn2", "Elavl4", "Actl6b", "Zcchc12", "Aplp1", "Stmn4", "Rtn1"),
  BranchA = c("Etv1", "Tbx3", "Pbx3", "Nos1", "Gal", "Ntng1", "Calb1"),
  BranchB = c("Bnc2", "Ndufa4l2", "Dlx5", "Cck", "Nmu")
)

save_featureplots_to_pdf(
  E16_barcodes_E3, cell_type_features,
  output_pdf = file.path(fig_dir_E3, "E16_embryo3_featureplots.pdf"),
  reduction = "umap.harmony",
  pt.size = 0.5
)

save_featureplots_to_pdf(
  E16_barcodes_E3, ENS_features,
  output_pdf = file.path(fig_dir_E3, "E16_embryo3_ENSrelated_featureplots.pdf"),
  reduction = "umap.harmony",
  pt.size = 0.5
)

# Label transfer from Experiment 1 (requires annotated E1 from Expt1_Processing.Rmd)
E1 <- readRDS("Data/E16_barcodes_E1.rds")
if (!"assigned_cell_type" %in% colnames(E1@meta.data)) {
  stop(
    "E1 object lacks 'assigned_cell_type'. Run Expt1_Processing.Rmd first and save ",
    "the annotated object to Data/E16_barcodes_E1.rds."
  )
}

E16_barcodes_E3.anchors <- FindTransferAnchors(reference = E1, query = E16_barcodes_E3, dims = 1:30,reference.reduction = "pca", normalization.method = "SCT",recompute.residuals = FALSE)
E16_barcodes_E3_predictions <- TransferData(anchorset = E16_barcodes_E3.anchors, refdata = E1$assigned_cell_type, dims = 1:30)
E16_barcodes_E3 <- AddMetaData(E16_barcodes_E3, metadata = E16_barcodes_E3_predictions)
table(E16_barcodes_E3_predictions$predicted.id)
DimPlot(E16_barcodes_E3, reduction = "umap.harmony", group.by = 'predicted.id', label=T,repel=T)
save_fig_E3("predictedid_E3.pdf", width = 10, height = 6)

E16_barcodes_E3$predicted.clusters <- ifelse(E16_barcodes_E3$prediction.score.max >= 0.5, 
                                             E16_barcodes_E3$predicted.id, NA)
DimPlot(E16_barcodes_E3, label = TRUE, repel = TRUE, reduction = "umap.harmony", group.by = 'predicted.clusters',split.by = 'sample')
save_fig_E3("predictedid_E3_with50threshold.pdf", width = 10, height = 6)


E16_barcodes_E3$Tomato <- GetAssayData(E16_barcodes_E3, layer = "counts")["Tomato-N", ] > 0

E16_barcodes_E3$Tomato <- factor(
  E16_barcodes_E3$Tomato,
  levels = c(FALSE, TRUE),
  labels = c("Tomato−", "Tomato+")
)
table(E16_barcodes_E3$Tomato) # sSupplementary 16752 /18023 (92.9 %unique()z%)
FeaturePlot(E16_barcodes_E3, features = c("Tomato-N"), reduction = "umap.harmony", order = T) + xlab("UMAP1") + ylab("UMAP2")

E16_barcodes_E3 <- PrepSCTFindMarkers(E16_barcodes_E3)
E3_markers <- FindAllMarkers(E16_barcodes_E3, only.pos = TRUE, test.use = "wilcox", min.pct = 0.25, logfc.threshold = 0.25, group.by = "harmony_clusters")

top50_E3_markers <- E3_markers %>%
  group_by(cluster) %>%
  top_n(50, avg_log2FC) %>%
  ungroup() %>%
  select(cluster, everything())

write.csv(
  top50_E3_markers,
  file = "Data/top50_E3_markers.csv",
  row.names = FALSE
)

E16_barcodes_E3 <- subset(E16_barcodes_E3, idents = "5", invert = TRUE)
DimPlot(E16_barcodes_E3, reduction = "umap.harmony", label = TRUE)

DefaultAssay(E16_barcodes_E3)='RNA'
# optional: remove old normalized/scaled layers
E16_barcodes_E3[["RNA"]]$data.E16_stomach_barcodes_E3<- NULL
E16_barcodes_E3[["RNA"]]$scale.data.E16_stomach_barcodes_E3<- NULL
E16_barcodes_E3[["RNA"]]$data.E16_jejileum_barcodes_E3<- NULL
E16_barcodes_E3[["RNA"]]$scale.data.E16_jejileum_barcodes_E3<- NULL
E16_barcodes_E3[["RNA"]]$data.E16_colon_barcodes_E3<- NULL
E16_barcodes_E3[["RNA"]]$scale.data.E16_colon_barcodes_E3<- NULL

E16_barcodes_E3[["RNA"]] <- split(E16_barcodes_E3[["RNA"]], f = E16_barcodes_E3$sample)
E16_barcodes_E3[["SCT"]] <- NULL

E16_barcodes_E3 <- SCTransform(E16_barcodes_E3, vars.to.regress = c("percent.mt", "nGene","nUMI"), verbose = TRUE)
VarGenes_E16_barcodes_E3 <- VariableFeatures(object = E16_barcodes_E3)
VarGenes_E16_barcodes_E3 <- VarGenes_E16_barcodes_E3[!VarGenes_E16_barcodes_E3 %in% c("Xist", "Gm13305", "Tsix", "Gm8730", "Eif2s3y", "Ddx3y", "Uty", "Kdm5d")]
E16_barcodes_E3 <- RunPCA(E16_barcodes_E3, verbose = FALSE, features = VarGenes_E16_barcodes_E3)
E16_barcodes_E3 <- RunUMAP(E16_barcodes_E3, reduction = "pca", dims = 1:50, n.neighbors = 55L, min.dist = 0.40, n.epochs = 1000, local.connectivity = 1L, seed.use = 111)
DimPlot(E16_barcodes_E3, reduction = "umap", group.by = "sample")

# Integration using Harmony
E16_barcodes_E3[["harmony"]] <- NULL
E16_barcodes_E3[["umap.harmony"]] <- NULL
E16_barcodes_E3 <- IntegrateLayers(
  object = E16_barcodes_E3, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  assay = "SCT", verbose = FALSE, group.by.vars = "sample"
)
ElbowPlot_harmony(E16_barcodes_E3, ndims = 50, reduction = "harmony") # theta = 3
E16_barcodes_E3 <- RunUMAP(E16_barcodes_E3, reduction = "harmony", dims = 1:50, reduction.name = "umap.harmony", random.seed = 222,n.neighbors = 50L, min.dist = 0.10, n.epochs = 1000) # set neighbour =30, SCP disppeared
E16_barcodes_E3 <- FindNeighbors(E16_barcodes_E3, reduction = "harmony", dims = 1:50) # set neighbour =30, SCP disappeared
E16_barcodes_E3 <- FindClusters(E16_barcodes_E3, resolution = 2.2, cluster.name = "harmony_clusters", algorithm = 4, random.seed = 222)
DimPlot(E16_barcodes_E3, reduction = "umap.harmony", group.by = "sample")
DimPlot(E16_barcodes_E3, reduction = "umap.harmony", label = TRUE)


E16_barcodes_E3 <- subset(E16_barcodes_E3, idents = "16", invert = TRUE)
DimPlot(E16_barcodes_E3, reduction = "umap.harmony", label = TRUE)

DefaultAssay(E16_barcodes_E3)='RNA'
E16_barcodes_E3[["RNA"]] <- split(E16_barcodes_E3[["RNA"]], f = E16_barcodes_E3$sample)
E16_barcodes_E3[["SCT"]] <- NULL

E16_barcodes_E3 <- SCTransform(E16_barcodes_E3, vars.to.regress = c("percent.mt", "nGene","nUMI"), verbose = TRUE)
VarGenes_E16_barcodes_E3 <- VariableFeatures(object = E16_barcodes_E3)
VarGenes_E16_barcodes_E3 <- VarGenes_E16_barcodes_E3[!VarGenes_E16_barcodes_E3 %in% c("Xist", "Gm13305", "Tsix", "Gm8730", "Eif2s3y", "Ddx3y", "Uty", "Kdm5d")]
E16_barcodes_E3 <- RunPCA(E16_barcodes_E3, verbose = FALSE, features = VarGenes_E16_barcodes_E3)
E16_barcodes_E3 <- RunUMAP(E16_barcodes_E3, reduction = "pca", dims = 1:50, n.neighbors = 55L, min.dist = 0.40, n.epochs = 1000, local.connectivity = 1L, seed.use = 111)
DimPlot(E16_barcodes_E3, reduction = "umap", group.by = "sample")

# Integration using Harmony
E16_barcodes_E3[["harmony"]] <- NULL
E16_barcodes_E3[["umap.harmony"]] <- NULL
E16_barcodes_E3 <- IntegrateLayers(
  object = E16_barcodes_E3, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  assay = "SCT", verbose = FALSE, group.by.vars = "sample"
)
ElbowPlot_harmony(E16_barcodes_E3, ndims = 50, reduction = "harmony") # theta = 3
E16_barcodes_E3 <- RunUMAP(E16_barcodes_E3, reduction = "harmony", dims = 1:50, reduction.name = "umap.harmony",n.neighbors = 60L, min.dist = 0.1, n.epochs = 500) # set neighbour =30, SCP disppeared
E16_barcodes_E3 <- FindNeighbors(E16_barcodes_E3, reduction = "harmony", dims = 1:50) # set neighbour =30, SCP disappeared
E16_barcodes_E3 <- FindClusters(E16_barcodes_E3, resolution = 1.2, cluster.name = "harmony_clusters", algorithm = 4, random.seed = 222)
DimPlot(E16_barcodes_E3, reduction = "umap.harmony", group.by = "sample")
DimPlot(E16_barcodes_E3, reduction = "umap.harmony", label = TRUE)


# remove cluster 
E16_barcodes_E3 <- subset(E16_barcodes_E3, idents = "9", invert = TRUE)
DimPlot(E16_barcodes_E3, reduction = "umap.harmony", label = TRUE)

DefaultAssay(E16_barcodes_E3)='RNA'

E16_barcodes_E3[["SCT"]] <- NULL

E16_barcodes_E3 <- SCTransform(E16_barcodes_E3, vars.to.regress = c("percent.mt", "nGene","nUMI"), verbose = TRUE)
VarGenes_E16_barcodes_E3 <- VariableFeatures(object = E16_barcodes_E3)
VarGenes_E16_barcodes_E3 <- VarGenes_E16_barcodes_E3[!VarGenes_E16_barcodes_E3 %in% c("Xist", "Gm13305", "Tsix", "Gm8730", "Eif2s3y", "Ddx3y", "Uty", "Kdm5d")]
E16_barcodes_E3 <- RunPCA(E16_barcodes_E3, verbose = FALSE, features = VarGenes_E16_barcodes_E3)
E16_barcodes_E3 <- RunUMAP(E16_barcodes_E3, reduction = "pca", dims = 1:50, n.neighbors = 55L, min.dist = 0.40, n.epochs = 1000, local.connectivity = 1L, seed.use = 111)
DimPlot(E16_barcodes_E3, reduction = "umap", group.by = "sample")

# Integration using Harmony
E16_barcodes_E3[["harmony"]] <- NULL
E16_barcodes_E3[["umap.harmony"]] <- NULL
E16_barcodes_E3 <- IntegrateLayers(
  object = E16_barcodes_E3, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  assay = "SCT", verbose = FALSE, group.by.vars = "sample"
)
ElbowPlot_harmony(E16_barcodes_E3, ndims = 50, reduction = "harmony") # theta = 3
E16_barcodes_E3 <- RunUMAP(E16_barcodes_E3, reduction = "harmony", dims = 1:50, reduction.name = "umap.harmony",n.neighbors = 50L, min.dist = 0.1, n.epochs = 1000) # set neighbour =30, SCP disppeared
E16_barcodes_E3 <- FindNeighbors(E16_barcodes_E3, reduction = "harmony", dims = 1:50) # set neighbour =30, SCP disappeared
E16_barcodes_E3 <- FindClusters(E16_barcodes_E3, resolution = 1.1, cluster.name = "harmony_clusters", algorithm = 4, random.seed = 222)
DimPlot(E16_barcodes_E3, reduction = "umap.harmony", group.by = "sample")
DimPlot(E16_barcodes_E3, reduction = "umap.harmony", label = TRUE)

saveRDS(object = E16_barcodes_E3, "Data/E16_barcodes_E3.rds")

# Exporting cellIDs for TREX analysis
ids_E3 <- rownames(E16_barcodes_E3@meta.data)
head(ids_E3)
ids_E3 <- gsub("-1$", "", ids_E3)
head(ids_E3)
indexed_ids_E3 <- data.frame(
  index = seq(0, length(ids_E3) - 1),
  cell_id = ids_E3
)

write.table(
  indexed_ids_E3,
  file = "Data/cellids_E3.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)











