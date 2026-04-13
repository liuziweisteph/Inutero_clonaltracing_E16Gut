library(Seurat)
library(sctransform)
library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(SoupX)
library(DropletUtils)
library(scCustomize)
library(scran)
source("Helpers/SoupXDoubletfinder.R")

# Comparing background ambient RNA Stating marker genes for diverse cell types in the developing mouse gut
marker_genes <- c(
  "Myh11", "Acta2", "Actg2", "Col1a1", "Lum", "Cdh1", "Epcam", "Cd36", "Pecam1", "Upk3b", "Msln", "Wt1",
  "Phox2b", "Hand2", "Ret", "Kit", "Ano1", "Lyz2", "Mrc1", "Trbc2", "Cd52", "Ptprc", "Rgs5", "Reln", "Pdgfrb", "Pdgfra"
)

external_data_dir <- "~/scRNA_seq/Trex_ENS/Submission/Data"  # Edit this path to your local raw data location
in_data_dir <- file.path(external_data_dir, "counts")
in_data_folders <- list.dirs(in_data_dir, recursive = TRUE, full.names = TRUE)
raw_data_folders <- in_data_folders[grepl("raw_feature_bc_matrix$", in_data_folders)]

# Match raw and filtered by parent sample folder
sample_names <- basename(dirname(raw_data_folders))

pairs <- data.frame(sample = sample_names, raw = raw_data_folders, filtered = sub(
  "raw_feature_bc_matrix", "filtered_feature_bc_matrix",
  raw_data_folders
), stringsAsFactors = FALSE)

pairs

ambient_summary <- analyze_ambient(pairs, marker_genes)
write.csv(ambient_summary, "Data/ambient_RNA_QC_summary.csv", row.names = FALSE)

# After running SoupX
in_data_dir_afterSoupX <- file.path(external_data_dir, "SoupX")
in_data_folders_afterSoupX <- list.dirs(in_data_dir_afterSoupX, recursive = TRUE, full.names = TRUE)
raw_data_folders_afterSoupX <- in_data_dir_afterSoupX[grepl("raw_feature_bc_matrix$", in_data_folders_afterSoupX)]

# Match raw and filtered by parent sample folder
sample_names_afterSoupX <- basename(dirname(raw_data_folders_afterSoupX))

pairs_afterSoupX <- data.frame(sample = sample_names_afterSoupX, raw = raw_data_folders_afterSoupX, filtered = sub(
  "raw_feature_bc_matrix",
  "filtered_feature_bc_matrix", raw_data_folders_afterSoupX
), stringsAsFactors = FALSE)

pairs_afterSoupX

ambient_summary_afterSoupX <- analyze_ambient(pairs_afterSoupX, marker_genes)
write.csv(ambient_summary_afterSoupX, "Data/ambient_RNA_QC_summary_afterSoupX.csv",
  row.names = FALSE
)
