# Author: Ziwei Liu, Karolinska Institutet
#
# SoupX-based functions below are adapted from the SoupX package and workflow described in:
# Matthew D Young, Sam Behjati. SoupX removes ambient RNA contamination from droplet-based
# single-cell RNA sequencing data (Open Access). GigaScience, Volume 9, Issue 12, December 2020,
# giaa151; bioRxiv 303727. https://doi.org/10.1093/gigascience/giaa151
# https://academic.oup.com/gigascience/article/9/12/giaa151/6049831
#
# DoubletFinder-based functions are adapted from the DoubletFinder package and:
# McGinnis, C.S., Murrow, L.M., Gartner, Z.J. (2019). DoubletFinder: Doublet Detection in
# Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors. Cell Systems 8(4), 329-337.
# https://doi.org/10.1016/j.cels.2019.03.003

# --- SoupX (Young & Behjati 2020) ---
analyze_ambient <- function(pairs, marker_genes) {
  ambient_results <- list()

  for (i in seq_len(nrow(pairs))) {
    dataset_name <- pairs$sample[i]
    raw_path <- pairs$raw[i]
    filtered_path <- pairs$filtered[i]

    message("Processing: ", dataset_name)

    tod <- Read10X(raw_path)
    toc <- Read10X(filtered_path)
    soup <- SoupChannel(tod, toc)
    empty <- tod[, colSums(tod) < 100, drop = FALSE] # Empty droplets
    bg_profile <- rowMeans(empty) # background
    # Marker leakage into empty droplets
    markers_present <- intersect(marker_genes, rownames(empty))
    marker_bg <- rowMeans(empty[markers_present, , drop = FALSE])
    mean_marker_bg <- mean(marker_bg)

    # Estimate contamination fraction
    clusters <- quickCluster(soup$toc)
    soup <- setClusters(soup, clusters)
    soup <- autoEstCont(soup)
    contam_frac <- soup$metaData$rho

    # Summaries
    result <- data.frame(
      Dataset = dataset_name,
      N_droplets = ncol(tod),
      N_cells = ncol(toc),
      Mean_background_UMI = mean(bg_profile),
      Mean_marker_background = mean_marker_bg,
      Mean_contamination_fraction = mean(contam_frac),
      Median_contamination_fraction = median(contam_frac),
      stringsAsFactors = FALSE
    )

    # Store per-dataset outputs
    ambient_results[[dataset_name]] <- list(
      result = result, bg_profile = bg_profile, contam_frac = contam_frac
    )
  }

  # Combine summary table
  ambient_summary <- do.call(
    rbind, lapply(ambient_results, `[[`, "result")
  )

  return(ambient_summary)
}

# SoupX: correct counts and write 10x output (Young & Behjati 2020)
soupX_correction_counts_andsave <- function(in_data_dir, out_parent_dir) {
  in_data_folders <- list.dirs(in_data_dir, full.names = TRUE, recursive = FALSE)
  sample_dirs <- in_data_folders[
    sapply(in_data_folders, function(x) {
      file.exists(file.path(x, "raw_feature_bc_matrix")) &&
        file.exists(file.path(x, "filtered_feature_bc_matrix"))
    })
  ]

  for (sample_path in sample_dirs) {
    sample_id <- basename(sample_path)
    message(paste0("Processing ", sample_id))
    raw_path <- file.path(sample_path, "raw_feature_bc_matrix")
    filtered_path <- file.path(sample_path, "filtered_feature_bc_matrix")

    sobj_path <- Read10X(data.dir = filtered_path)
    sobj <- CreateSeuratObject(counts = sobj_path, project = sample_id)
    sobj <- NormalizeData(sobj, verbose = FALSE)
    sobj <- FindVariableFeatures(object = sobj, nfeatures = 3000, verbose = FALSE, selection.method = "vst")
    sobj <- ScaleData(sobj, verbose = FALSE)
    sobj <- RunPCA(sobj, npcs = 30, verbose = FALSE)
    sobj <- FindNeighbors(sobj, dims = 1:30, verbose = FALSE)
    sobj <- FindClusters(sobj, resolution = 0.5, verbose = FALSE)
    sobj$soup_group <- sobj@meta.data[["seurat_clusters"]]

    raw <- Read10X(data.dir = raw_path)
    sc <- SoupChannel(raw, sobj[["RNA"]]$counts)
    sc <- setClusters(sc, sobj$soup_group)
    sc <- autoEstCont(sc, doPlot = FALSE)
    print(head(sc$soupProfile[order(sc$soupProfile$est, decreasing = T), ], n = 20))
    corrected_counts <- adjustCounts(sc)

    # Save
    output_dir <- file.path(out_parent_dir, paste0(sample_id, "_SoupX_corrected"))
    # dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    DropletUtils::write10xCounts(output_dir, corrected_counts, version = "3")
    message(paste("Saved corrected matrix for", sample_id, "to", output_dir))
  }
}


process_seurat_objects_fordb <- function(seu_list) {
  # Load Seurat objects
  for (i in seq_along(seu_list)) {
    seuobj <- seu_list[[i]]
    sample_id <- names(seu_list)[i]
    message(paste("Start loading from ", sample_id))

    # Normalisation and scaling
    seuobj <- NormalizeData(seuobj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    seuobj <- FindVariableFeatures(object = seuobj, nfeatures = 3000, verbose = FALSE, selection.method = "vst")
    seuobj <- ScaleData(seuobj, verbose = FALSE)
    seuobj <- RunPCA(seuobj, npcs = 50, verbose = FALSE)
    seuobj <- RunUMAP(seuobj, reduction = "pca", dims = 1:50, n.neighbors = 55L, min.dist = 0.40, n.epochs = 1000, local.connectivity = 1L)
    seuobj <- FindNeighbors(seuobj, dims = 1:50, verbose = FALSE)
    seuobj <- FindClusters(seuobj, resolution = 0.8, verbose = FALSE)
    seuobj$db_group <- Idents(seuobj)

    # Store in the list
    seu_list[[i]] <- seuobj
  }

  return(seu_list)
}

# --- DoubletFinder (McGinnis et al. 2019) ---
do_doubletfinder <- function(seu_list, expectRate) {
  for (i in seq_along(seu_list)) {
    seuobj <- seu_list[[i]]
    sample_id <- names(seu_list)[i]

    message(paste0("Performing Doubletfinder for ", sample_id))

    ## pK Identification (no ground-truth)
    sweep.res.list <- paramSweep(seuobj, PCs = 1:10, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    ## Optimal pK value
    optimal_pK <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric), "pK"]))

    ## Homotypic Doublet Proportion Estimate
    homotypic.prop <- modelHomotypic(seuobj$seurat_clusters)

    ## Number of expected doublets (adjust based on your expected rate)
    nExp_poi <- round(expectRate[i] * ncol(seuobj)) # Expected rate based on cell recovery target
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

    ## Run DoubletFinder with adjusted parameters
    seuobj <- doubletFinder(seuobj,
      PCs = 1:30,
      pN = 0.25,
      pK = optimal_pK,
      nExp = nExp_poi.adj,
      reuse.pANN = NULL,
      sct = TRUE
    )

    # Extract doublet classifications
    classification_col <- grep("DF.classifications", colnames(seuobj@meta.data), value = TRUE)
    seuobj$db_class <- seuobj@meta.data[[classification_col]]

    # Store the Seurat objects in a list
    seu_list[[i]] <- seuobj
    rm(sweep.res.list, sweep.stats, bcmvn, optimal_pK, homotypic.prop, nExp_poi, nExp_poi.adj)
  }

  # Return seurat list
  return(seu_list)
}

# DoubletFinder: visualise DF classifications (McGinnis et al. 2019)
plot_doublets <- function(seu_list, group.by = "db_class", fig_dir) {
  # Ensure the directory exists
  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir, recursive = TRUE)
    message("Created directory: ", fig_dir)
  }

  for (i in seq_along(seu_list)) {
    seuobj <- seu_list[[i]]
    sample_id <- names(seu_list)[i]
    output_pdf <- file.path(fig_dir, paste0(sample_id, "_doublets_plots.pdf"))

    message("Generating plots for sample: ", sample_id)

    # Create plots
    p1 <- DimPlot(seuobj,
      reduction = "umap", group.by = group.by, pt.size = 0.5,
      label = FALSE, label.size = 5
    )

    p2 <- VlnPlot(
      seuobj,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
      pt.size = 0.2,
      group.by = group.by
    ) + ggtitle(paste("Sample:", sample_id))

    # Save to PDF
    pdf(output_pdf, width = 8, height = 8)
    print(plot_grid(p1, p2, ncol = 1))
    dev.off()

    message("Plots saved to: ", output_pdf)
  }
}


ElbowPlot_harmony <- function(object, reduction = "harmony", ndims = 50) {
  # Extract embeddings
  harmony_embeddings <- Embeddings(object, reduction = reduction)

  # Limit to requested number of dimensions
  ndims <- min(ndims, ncol(harmony_embeddings))
  harmony_embeddings <- harmony_embeddings[, 1:ndims, drop = FALSE]

  # Compute variance explained
  var_per_comp <- apply(harmony_embeddings, 2, var)
  var_explained <- var_per_comp / sum(var_per_comp) * 100

  # Create dataframe for ggplot
  df <- data.frame(
    Dimension = 1:ndims,
    Variance = var_explained[1:ndims]
  )

  # ggplot elbow plot
  ggplot(df, aes(x = Dimension, y = Variance)) +
    geom_point(size = 2, color = "steelblue") +
    geom_line(color = "steelblue") +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Elbow Plot for", reduction),
      x = paste0(toupper(substring(reduction, 1, 1)), substring(reduction, 2), " Component"),
      y = "Percent Variance Explained"
    )
}

save_all_clusters_to_one_pdf <- function(
  seurat_obj, unique_cell_types,
  out_pdf = "All_Datasets_Cluster_Highlight_Plots.pdf",
  highlight_color = "gold",
  background_color = "lightgray",
  reduction = "umap",
  pt.size = 1
) {
  # Start multi-page PDF
  pdf(out_pdf, width = 6, height = 5)

  for (cell_type in unique_cell_types) {
    p <- Cluster_Highlight_Plot(
      seurat_object    = seurat_obj,
      cluster_name     = cell_type,
      highlight_color  = "gold",
      background_color = "lightgray",
      pt.size          = 1,
      reduction        = "umap.harmony"
    ) +
      ggplot2::ggtitle(paste("Cluster Highlight:", cell_type))

    print(p) # patchwork/ggplot requires print() for PDF output
  }

  dev.off()
  message("Combined PDF saved")
}

save_featureplots_to_pdf <- function(
  seurat_obj, # single Seurat object
  feature_list, # named list of feature sets per cell type
  output_pdf = "E16_all_featureplots.pdf",
  pt.size = 0.5,
  reduction = "umap"
) {
  pdf(output_pdf, width = 6, height = 5, onefile = TRUE)

  for (cell_type in names(feature_list)) {
    features <- feature_list[[cell_type]]

    for (feature in features) {
      p <- FeaturePlot(
        object = seurat_obj,
        features = feature,
        pt.size = pt.size,
        reduction = reduction,
        order = TRUE
      ) + ggtitle(paste(cell_type, "-", feature))

      print(p) # required for patchwork/ggplot objects to render in PDF
    }
  }

  dev.off()
  message("All feature plots saved in: ", output_pdf)
}
