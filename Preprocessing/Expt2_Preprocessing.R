# Preprocessing Experiment 2: SoupX correction, QC, doublet removal, merge, SCTransform, Harmony
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
source("Helpers/SoupXDoubletfinder.R")

# External data directory — 
external_data_dir <- "~/scRNA_seq/Trex_ENS/Submission/Data"

# SoupX -------------------------------------------------------------------
in_data_dir <- file.path(external_data_dir, "counts/E2")
in_data_folders <- list.dirs(in_data_dir, recursive = TRUE, full.names = TRUE)
filtered_data_folders <- in_data_folders[grepl("filtered_feature_bc_matrix$", in_data_folders)]
raw_data_folders <- in_data_folders[grepl("raw_feature_bc_matrix$", in_data_folders)]
out_dir <- file.path(external_data_dir, "SoupX/E2")

soupX_correction_counts_andsave(in_data_dir, out_dir)

# Experiment 2 ---------------------------------------------------------------

X10_e16_s1_E2.data <- Read10X(data.dir = file.path(external_data_dir, "SoupX/E2/E16_stomach_SoupX_corrected/filtered_feature_bc_matrix"))
e16_s1_E2 <- CreateSeuratObject(counts = X10_e16_s1_E2.data, min.cells = 3, min.features = 200, project = "E16_stomach_barcodes_E2")
e16_s1_E2 <- PercentageFeatureSet(e16_s1_E2, pattern = "^mt-", col.name = "percent.mt")
VlnPlot(object = e16_s1_E2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
e16_s1_E2 <- subset(e16_s1_E2, subset = nFeature_RNA > 2000 & nFeature_RNA < 7500 & percent.mt < 10 & nCount_RNA > 2000 & nCount_RNA < 40000) # 10553 cells
VlnPlot(object = e16_s1_E2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
rm(X10_e16_s1_E2.data)

X10_e16_s2_E2.data <- Read10X(data.dir = file.path(external_data_dir, "SoupX/E2/E16_jejile_SoupX_corrected/filtered_feature_bc_matrix"))
e16_s2_E2 <- CreateSeuratObject(counts = X10_e16_s2_E2.data, min.cells = 3, min.features = 200, project = "E16_jejileum_barcodes_E2")
e16_s2_E2 <- PercentageFeatureSet(e16_s2_E2, pattern = "^mt-", col.name = "percent.mt")
VlnPlot(object = e16_s2_E2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
e16_s2_E2 <- subset(e16_s2_E2, subset = nFeature_RNA > 2000 & nFeature_RNA < 8000 & percent.mt < 10 & nCount_RNA > 2000 & nCount_RNA < 50000) # 9314 cells
VlnPlot(object = e16_s2_E2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
rm(X10_e16_s2_E2.data)

X10_e16_s3_E2.data <- Read10X(data.dir = file.path(external_data_dir, "SoupX/E2/E16_colon_SoupX_corrected/filtered_feature_bc_matrix"))
e16_s3_E2 <- CreateSeuratObject(counts = X10_e16_s3_E2.data, min.cells = 3, min.features = 200, project = "E16_colon_barcodes_E2")
e16_s3_E2 <- PercentageFeatureSet(e16_s3_E2, pattern = "^mt-", col.name = "percent.mt")
VlnPlot(object = e16_s3_E2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
e16_s3_E2 <- subset(e16_s3_E2, subset = nFeature_RNA > 2000 & nFeature_RNA < 7500 & percent.mt < 10 & nCount_RNA > 2000 & nCount_RNA < 50000) # 4002 cells
VlnPlot(object = e16_s3_E2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
rm(X10_e16_s3_E2.data)

E16_E2_list <- list("Stomach_E2" = e16_s1_E2, "Jej_ile_E2" = e16_s2_E2, "colon_E2" = e16_s3_E2)

# Doublet detection per sample ------------------------------------------------------
expectRate_E2 <- c(0.08, 0.072, 0.032)
names(expectRate_E2) <- c("Stomach_E2", "Jej_ile_E2", "colon_E2")

E16_E2_list <- process_seurat_objects_fordb(E16_E2_list)
E16_E2_list <- do_doubletfinder(E16_E2_list, expectRate_E2)

# Remove doublets ---------------------------------------------------------
E16_E2_list_clean <- sapply(E16_E2_list, function(seurat_obj) {
  subset(seurat_obj, subset = db_class == "Singlet")
}, simplify = FALSE)

cell_counts_E2 <- sapply(E16_E2_list_clean, function(seurat_obj) {
  ncol(seurat_obj)
})

print(cell_counts_E2)
# After SoupX: Stomach_E2   Jej_ile_E2    colon_E2
# 9760        8676        3884

# Merge datasets  ---------------------------------------------------------
E16_barcodes_E2 <- merge(E16_E2_list_clean[[1]], y = c(E16_E2_list_clean[[2]], E16_E2_list_clean[[3]]), add.cell.ids = c("e16_s1", "e16_ji1", "e16_c1"), project = "ENS_E16_barcodes_E2")

E16_barcodes_metadata_E2 <- E16_barcodes_E2@meta.data
# Add cell IDs to metadata
E16_barcodes_metadata_E2$cells <- rownames(E16_barcodes_metadata_E2)

# Rename columns
E16_barcodes_metadata_E2 <- E16_barcodes_metadata_E2 %>%
  dplyr::rename(
    seq_folder = orig.ident,
    nUMI = nCount_RNA,
    nGene = nFeature_RNA
  )

# Create sample column
E16_barcodes_metadata_E2$sample <- NA
E16_barcodes_metadata_E2$sample[which(str_detect(E16_barcodes_metadata_E2$cells, "^e16_s1_"))] <- "E16_stomach"
E16_barcodes_metadata_E2$sample[which(str_detect(E16_barcodes_metadata_E2$cells, "^e16_ji1_"))] <- "E16_jej_ileum"
E16_barcodes_metadata_E2$sample[which(str_detect(E16_barcodes_metadata_E2$cells, "^e16_c1_"))] <- "E16_colon"
E16_barcodes_metadata_E2$mitoRatio <- E16_barcodes_metadata_E2$percent.mt / 100
E16_barcodes_metadata_E2$log10GenesPerUMI <- log10(E16_barcodes_metadata_E2$nGene) / log10(E16_barcodes_metadata_E2$nUMI)
E16_barcodes_metadata_E2$batch <- "Embryo 2"
E16_barcodes_E2@meta.data <- E16_barcodes_metadata_E2

E16_barcodes_E2 <- SCTransform(E16_barcodes_E2, vars.to.regress = "percent.mt", verbose = TRUE)
VarGenes_E16_barcodes_E2 <- VariableFeatures(object = E16_barcodes_E2)
VarGenes_E16_barcodes_E2 <- VarGenes_E16_barcodes_E2[!VarGenes_E16_barcodes_E2 %in% c("Xist", "Gm13305", "Tsix", "Gm8730", "Eif2s3y", "Ddx3y", "Uty", "Kdm5d")]
E16_barcodes_E2 <- RunPCA(E16_barcodes_E2, verbose = FALSE, features = VarGenes_E16_barcodes_E2)
E16_barcodes_E2 <- RunUMAP(E16_barcodes_E2, reduction = "pca", dims = 1:50, n.neighbors = 55L, min.dist = 0.40, n.epochs = 1000, local.connectivity = 1L, seed.use = 111)
DimPlot(E16_barcodes_E2, reduction = "umap", group.by = "sample")

# Integration using Harmony
E16_barcodes_E2 <- IntegrateLayers(
  object = E16_barcodes_E2, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  assay = "SCT", verbose = FALSE, group.by.vars = "sample"
)
ElbowPlot_harmony(E16_barcodes_E2, ndims = 50, reduction = "harmony") # theta = 3
E16_barcodes_E2 <- RunUMAP(E16_barcodes_E2, reduction = "harmony", dims = 1:50, reduction.name = "umap.harmony", random.seed = 222) # With n.neighbors = 30, SCP cluster disappeared
E16_barcodes_E2 <- FindNeighbors(E16_barcodes_E2, reduction = "harmony", dims = 1:50) # With n.neighbors = 30, SCP cluster disappeared
E16_barcodes_E2 <- FindClusters(E16_barcodes_E2, resolution = 0.1, cluster.name = "harmony_clusters", algorithm = 4, random.seed = 222)
DimPlot(E16_barcodes_E2, reduction = "umap.harmony", group.by = "sample")
DimPlot(E16_barcodes_E2, reduction = "umap.harmony", label = TRUE)

VlnPlot(E16_barcodes_E2, features = c("nUMI", "nGene", "percent.mt"), ncol = 3, pt.size = 0.01)

genes <- c("Sox10", "Plp1", "Ascl1", "Elavl4", "Etv1", "Ndufa4l2", "Gfra3", "Meis2", "Lum", "Acta2")
FeaturePlot(E16_barcodes_E2, features = genes, order = T, reduction = "umap.harmony")
FeaturePlot(E16_barcodes_E2, features = c("Phox2b"), order = T, reduction = "umap.harmony")

# Reload cell cycle genes
mmus_s <- read.table("Data/mmus_s_genes.txt",
  sep = "\t", header = FALSE, stringsAsFactors = FALSE
)$V1
mmus_g2m <- read.table("Data/mmus_g2m_genes.txt",
  sep = "\t", header = FALSE, stringsAsFactors = FALSE
)$V1
mmus_s <- mmus_s[-1]
mmus_g2m <- mmus_g2m[-1]

E16_barcodes_E2 <- CellCycleScoring(E16_barcodes_E2,
  s.features = mmus_s,
  g2m.features = mmus_g2m, set.ident = TRUE
)

DimPlot(E16_barcodes_E2, group.by = "Phase", reduction = "umap.harmony")

Idents(E16_barcodes_E2) <- E16_barcodes_E2$harmony_clusters

cell_type_features <- list(
  Muscle = c("Myh11", "Acta2", "Actg2"),
  Fibroblasts = c("Col1a1", "Lum"),
  Epithelial = c("Cdh1", "Epcam"),
  Endothelial = c("Cd36", "Pecam1"),
  Mesothelial = c("Upk3b", "Msln", "Wt1"),
  ENS = c("Phox2b", "Hand2", "Ret"),
  ICC = c("Kit", "Ano1"),
  Immune = c("Lyz2", "Mrc1", "Trbc2", "Cd52", "Ptprc"),
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
  E16_barcodes_E2, cell_type_features,
  output_pdf = "Preprocessing/E16_embyro2_all_featureplots_finalised.pdf",
  reduction = "umap.harmony",
  pt.size = 0.5
)

save_featureplots_to_pdf(
  E16_barcodes_E2, ENS_features,
  output_pdf = "Preprocessing/E16_embyro2_ENSrelated_featureplots_finalised.pdf",
  reduction = "umap.harmony",
  pt.size = 0.5
)

saveRDS(object = E16_barcodes_E2, "Data/E16_barcodes_E2.rds")
