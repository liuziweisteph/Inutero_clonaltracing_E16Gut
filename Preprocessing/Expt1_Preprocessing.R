# Preprocessing Experiment 1: QC, doublet removal, merge, SCTransform, Harmony integration
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
library(scCustomize)
library(biomaRt)
py_module_available(module = "leidenalg")
source("Helpers/SoupXDoubletfinder.R")

# External data directory — 
external_data_dir <- "~/scRNA_seq/Trex_ENS/Submission/Data"

# Experiment 1 ---------------------------------------------------------------

X10_e16_s1.data <- Read10X(data.dir = file.path(external_data_dir, "counts/E1/E16_stomach_E1"))
e16_s1 <- CreateSeuratObject(counts = X10_e16_s1.data, min.cells = 3, min.features = 200, project = "E16_stomach_barcodes")
rm(X10_e16_s1.data)
e16_s1 <- PercentageFeatureSet(e16_s1, pattern = "^mt-", col.name = "percent.mt")
VlnPlot(object = e16_s1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
e16_s1 <- subset(e16_s1, subset = nFeature_RNA > 2000 & nFeature_RNA < 9000 & percent.mt < 10 & nCount_RNA > 2000 & nCount_RNA < 50000) # 4142 cells
VlnPlot(object = e16_s1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)

X10_e16_s2.data <- Read10X(data.dir = file.path(external_data_dir, "counts/E1/E16_ile_E1"))
e16_s2 <- CreateSeuratObject(counts = X10_e16_s2.data, min.cells = 3, min.features = 200, project = "E16_ileum_barcodes")
rm(X10_e16_s2.data)
e16_s2 <- PercentageFeatureSet(e16_s2, pattern = "^mt-", col.name = "percent.mt")
VlnPlot(object = e16_s2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
e16_s2 <- subset(e16_s2, subset = nFeature_RNA > 2000 & nFeature_RNA < 9000 & percent.mt < 10 & nCount_RNA > 2000 & nCount_RNA < 50000) # 1426 cells
VlnPlot(object = e16_s2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)

E16_E1_list <- list("Stomach" = e16_s1, "Ileum" = e16_s2)

# Doublet detection per sample
expectRate <- c(0.032, 0.016)
names(expectRate) <- c("Stomach", "Ileum")

E16_E1_list <- process_seurat_objects_fordb(E16_E1_list)
E16_E1_list <- do_doubletfinder(E16_E1_list, expectRate)

# Visualizing doublets identified by DoubletFinder
plot_doublets(E16_E1_list, fig_dir = "Preprocessing/Doubletfinder/")

# Removing the doublets from the dataset
E16_E1_list_clean <- sapply(E16_E1_list, function(seurat_obj) {
  subset(seurat_obj, subset = db_class == "Singlet")
}, simplify = FALSE)

# Merge datasets  ---------------------------------------------------------
E16_barcodes <- merge(E16_E1_list_clean[[1]], y = c(E16_E1_list_clean[[2]]), add.cell.ids = c("e16_s1", "e16_i1"), project = "ENS_E16_barcodes")
E16_barcodes_metadata <- E16_barcodes@meta.data

# Add cell IDs to metadata
E16_barcodes_metadata$cells <- rownames(E16_barcodes_metadata)

# Rename columns
E16_barcodes_metadata <- E16_barcodes_metadata %>%
  dplyr::rename(
    seq_folder = orig.ident,
    nUMI = nCount_RNA,
    nGene = nFeature_RNA
  )

# Create sample column
E16_barcodes_metadata$sample <- NA
E16_barcodes_metadata$sample[which(str_detect(E16_barcodes_metadata$cells, "^e16_s1_"))] <- "E16_stomach"
E16_barcodes_metadata$sample[which(str_detect(E16_barcodes_metadata$cells, "^e16_i1_"))] <- "E16_ileum"

E16_barcodes_metadata$mitoRatio <- E16_barcodes_metadata$percent.mt / 100
E16_barcodes_metadata$log10GenesPerUMI <- log10(E16_barcodes_metadata$nGene) / log10(E16_barcodes_metadata$nUMI)

# Visualize the number of cell counts per sample
E16_barcodes_metadata %>%
  ggplot(aes(x = sample, fill = sample)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("E16_barcodes ")

E16_barcodes@meta.data <- E16_barcodes_metadata
E16_barcodes <- SCTransform(E16_barcodes, vars.to.regress = c("percent.mt"), verbose = TRUE)
VarGenes_E16_barcodes <- VariableFeatures(object = E16_barcodes)
VarGenes_E16_barcodes <- VarGenes_E16_barcodes[!VarGenes_E16_barcodes %in% c("Xist", "Gm13305", "Tsix", "Gm8730", "Eif2s3y", "Ddx3y", "Uty", "Kdm5d")]
E16_barcodes <- RunPCA(E16_barcodes, verbose = FALSE, features = VarGenes_E16_barcodes)
E16_barcodes <- RunUMAP(E16_barcodes, reduction = "pca", dims = 1:50, n.neighbors = 55L, min.dist = 0.40, n.epochs = 1000, local.connectivity = 1L, seed.use = 111)
DimPlot(E16_barcodes, reduction = "umap", group.by = "sample")

# Integration using Harmony
E16_barcodes <- IntegrateLayers(
  object = E16_barcodes, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  assay = "SCT", verbose = FALSE
)
E16_barcodes <- RunUMAP(E16_barcodes, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", random.seed = 11)
E16_barcodes <- FindNeighbors(E16_barcodes, reduction = "harmony", dims = 1:30)
E16_barcodes <- FindClusters(E16_barcodes, resolution = 0.05, cluster.name = "harmony_clusters", algorithm = 4, random.seed = 111)
DimPlot(E16_barcodes, reduction = "umap.harmony", group.by = "sample")
DimPlot(E16_barcodes, reduction = "umap.harmony", label = TRUE)

E16_barcodes <- FindSubCluster(
  E16_barcodes,
  cluster = 1,
  graph.name = "SCT_snn",
  subcluster.name = "subcluster1", resolution = 0.4
)
DimPlot(E16_barcodes, reduction = "umap.harmony", label = TRUE, group.by = "subcluster1")

E16_barcodes <- PrepSCTFindMarkers(E16_barcodes)
subcluster_markers <- FindAllMarkers(E16_barcodes, only.pos = TRUE, test.use = "wilcox", min.pct = 0.25, logfc.threshold = 0.25, group.by = "subcluster")

top10_subcluster_markers <- subcluster_markers %>%
  group_by(cluster) %>%
  top_n(20, avg_log2FC)

# Decided to remove cluster 1_9 due to high mitochondrial content
FeaturePlot(E16_barcodes, features = c("Adam12", "Prss35", "Dpt", "Top2a", "Ebf1", "Ebf2", "Sox10", "Phox2b"), order = T, pt.size = 0.2, reduction = "umap.harmony")

E16_barcodes <- FindSubCluster(
  E16_barcodes,
  cluster = 2,
  graph.name = "SCT_snn",
  subcluster.name = "subcluster2",
  resolution = 0.2
)

DimPlot(E16_barcodes, reduction = "umap.harmony", label = TRUE, group.by = "subcluster2")

E16_barcodes$subcluster_combined <- NA

# Fill in subcluster1 info
if ("subcluster1" %in% colnames(E16_barcodes@meta.data)) {
  E16_barcodes$subcluster_combined[!is.na(E16_barcodes$subcluster1)] <-
    E16_barcodes$subcluster1[!is.na(E16_barcodes$subcluster1)]
}

# Identify cluster 2 cells
cluster2_cells <- Idents(E16_barcodes) == "2"

# Update only these cells with subcluster2 info
E16_barcodes$subcluster_combined[cluster2_cells] <- E16_barcodes$subcluster2[cluster2_cells]

# For cells with no subcluster info, keep their original cluster label
E16_barcodes$subcluster_combined[is.na(E16_barcodes$subcluster_combined)] <- as.character(E16_barcodes$seurat_clusters)

DimPlot(E16_barcodes, reduction = "umap.harmony", label = TRUE, group.by = "subcluster_combined")

Idents(E16_barcodes) <- E16_barcodes$subcluster_combined

# Decided to remove cluster 1_9 due to high mitochondrial content
E16_barcodes <- subset(E16_barcodes, idents = "1_9", invert = TRUE)
DimPlot(E16_barcodes, reduction = "umap.harmony", label = TRUE, group.by = "subcluster_combined")

# Cell cycle phase assignment
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

ensembl <- useMart("ensembl")
# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x) {
  require("biomaRt")
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  genesV2 <- getLDS(
    attributes = c("hgnc_symbol"),
    filters = "hgnc_symbol",
    values = x,
    mart = human,
    attributesL = c("mgi_symbol"),
    martL = mouse, uniqueRows = T
  )
  mousex <- unique(genesV2[, 2])
  human_genes_number <- length(x)
  mouse_genes_number <- length(mousex)
  if (human_genes_number != mouse_genes_number) {
    genes_not_trans <- setdiff(x, genesV2$HGNC.symbol)
    print("These genes could not be translated:")
    print(genes_not_trans)
    print(paste("A total number of ", length(genes_not_trans), "genes could not be translated!"), sep = " ")
  } else {
    print("All genes were translated successfully!")
  }
  return(mousex)
}

mmus_s <- convertHumanGeneList(s.genes) # "MLF1IP" "POLD3"  "ATAD2"  3 genes cannot be converted, 41 genes
mmus_g2m <- convertHumanGeneList(g2m.genes) # "TMPO"   "FAM64A" "HN1" 3 genes cannot be converted, 54 genes

E16_barcodes <- CellCycleScoring(E16_barcodes,
  s.features = mmus_s,
  g2m.features = mmus_g2m, set.ident = TRUE
)

write.table(mmus_s, file = "Data/mmus_s_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(mmus_g2m, file = "Data/mmus_g2m_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

DimPlot(E16_barcodes, group.by = "Phase", reduction = "umap.harmony")

Idents(E16_barcodes) <- E16_barcodes$harmony_clusters

cell_type_features <- list(
  Muscle = c("Myh11", "Acta2", "Actg2"),
  Fibroblasts = c("Col1a1", "Lum", "Vim"),
  Epithelial = c("Cdh1", "Epcam"),
  Endothelial = c("Cd36", "Pecam1"),
  Mesothelial = c("Upk3b", "Msln", "Wt1"),
  ENS = c("Phox2b", "Hand2", "Ret"),
  ICC = c("Kit", "Ano1"),
  Immune = c("Lyz2", "Mrc1", "Trbc2", "Cd52", "Ptprc"),
  Pericytes = c("Rgs5", "Reln", "Cox4i2", "Pdgfrb", "Pdgfra"), # Pdgfra negative
  Myofibroblast = c("Trpc3", "Acta2", "Tagln"),
  Proliferation = c("Top2a", "Mki67"),
  Enterocytes = c("Vil1", "Krt20"),
  Enteroendocrine = c("Chga", "Tph1"),
  Goblet = c("Fcgbp", "Clca1"),
  Stomach = c("Pitx1")
)

save_featureplots_to_pdf(
  E16_barcodes, cell_type_features,
  output_pdf = "Preprocessing/E16_expt1_all_featureplots.pdf",
  reduction = "umap.harmony",
  pt.size = 0.5
)

saveRDS(object = E16_barcodes, "Data/E16_barcodes_E1.rds")
