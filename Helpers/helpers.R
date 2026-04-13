# Function to create a list of Seurat objects from multiple input directories
create_seurat_objects <- function(in_data_folders, SoupX = TRUE, min.cells = 3, min.features = 200, project = "Developing ENS", fig_dir, verbose = TRUE) {
  seu_list <- list() # Initialize an empty list to store Seurat objects

  # Read 10X data from the SoupX_corrected folder or filtered_feature_bc_matrix folder
  data_folders <- if (SoupX) {
    in_data_folders[grepl("SoupX_corrected$", in_data_folders)]
  } else {
    in_data_folders[grepl("filtered_feature_bc_matrix$", in_data_folders)]
  }

  # Create output directory if it doesn't exist
  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir, recursive = TRUE)
    if (verbose) message(paste("Created figure directory:", fig_dir))
  }

  for (i in seq_along(data_folders)) {
    path <- data_folders[i]
    sample_id <- basename(dirname(path))

    # Print which sample is being processed
    print(paste0("Processing sample: ", sample_id))

    # Print the selected folder path
    print(paste0("Reading data from: ", path))

    # Read 10X data from the selected folder
    indata <- Read10X(data.dir = path)

    # Create Seurat object and compute mitochondrial percentage
    seuobj <- CreateSeuratObject(counts = indata, min.cells = min.cells, min.features = min.features, project = project) %>%
      PercentageFeatureSet(pattern = "^mt-", col.name = "percent.mt")

    # Add info to metadata
    seuobj$sample_id <- sample_id
    seuobj$time_point <- sub("_.*", "", sample_id)
    seuobj$region <- sub("^[^_]+_([^_]+)_.*$", "\\1", sample_id)

    output_pdf <- paste0(fig_dir, sample_id, "_QC_VlnPlots.pdf")
    p1 <- VlnPlot(object = seuobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01) + ggtitle(paste0(sample_id, " before filtering"))
    pdf(output_pdf, width = 8, height = 4)
    print(p1)
    dev.off()
    message(paste0(sample_id, " QC VlnPlots saved to: ", output_pdf))

    # Store in the list
    seu_list[[sample_id]] <- seuobj
  }

  print("All samples loaded successfully!")
  return(seu_list)
}

filter_seurat_objects <- function(seu_list,
                                  thresholds_list,
                                  fig_dir = "Figures/qc_figures",
                                  log_file = "Figures/qc_logs/filtering_summary.txt",
                                  save_individual = TRUE) {
  if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
  if (!dir.exists(dirname(log_file))) dir.create(dirname(log_file), recursive = TRUE)
  write("", file = log_file) # Clear log

  seu_list_new <- list()

  for (sample_id in names(seu_list)) {
    seuobj <- seu_list[[sample_id]]
    thresholds <- thresholds_list[[sample_id]]

    # Generate before-filtering violin plot
    Idents(seuobj) <- paste0(sample_id, " - Before Filtering")
    p_before <- VlnPlot(
      object = seuobj,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
      ncol = 3, pt.size = 0.01
    )

    # Record initial count
    n_before <- ncol(seuobj)

    # Apply filtering
    seuobj <- subset(seuobj,
      subset =
        nFeature_RNA > thresholds$min_genes & nFeature_RNA < thresholds$max_genes &
          nCount_RNA > thresholds$min_counts & nCount_RNA < thresholds$max_counts &
          percent.mt < thresholds$max_mt
    )

    n_after <- ncol(seuobj)

    # Generate after-filtering violin plot
    Idents(seuobj) <- paste0(sample_id, " - After Filtering")
    p_after <- VlnPlot(
      object = seuobj,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
      ncol = 3, pt.size = 0.01
    )

    # Combine and save to one PDF
    combined_plot <- p_before / p_after # patchwork layout (vertical)
    pdf_path <- file.path(fig_dir, paste0(sample_id, "_QC_VlnPlots.pdf"))
    ggsave(filename = pdf_path, plot = combined_plot, width = 10, height = 12)
    message(paste0("Violin plots saved to: ", pdf_path))

    # Log filtering summary
    filtered_out <- n_before - n_after
    pct_filtered <- round(100 * filtered_out / n_before, 1)
    msg <- paste0(
      "Sample: ", sample_id,
      " — Filtered out ", filtered_out, " cells (", pct_filtered, "%) ",
      "[nFeature_RNA: ", thresholds$min_genes, "-", thresholds$max_genes,
      ", nCount_RNA: ", thresholds$min_counts, "-", thresholds$max_counts,
      ", percent.mt < ", thresholds$max_mt, "]"
    )

    message(msg)
    cat(msg, "\n", file = log_file, append = TRUE)

    # Save result
    seu_list_new[[sample_id]] <- seuobj
    if (save_individual) assign(sample_id, seuobj, envir = .GlobalEnv)
  }

  return(seu_list_new)
}


process_seurat_objects_LogNormalize <- function(seu_list,
                                                var_genes_to_exclude = c("Xist", "Tsix", "Ddx3y", "Uty", "Eif2s3y", "Kdm5d", "Gm29650", "Malat1", "Gm42418", "AY036118", "Gm13305"),
                                                pca_dims = 50, umap_neighbors = 55L, umap_min_dist = 0.40, umap_epochs = 1000, umap_local_connectivity = 1L,
                                                seed = 11, cluster_resolution = 1.8, cluster_algorithm = 4, verbose = TRUE) {
  # Load Seurat objects
  for (i in seq_along(seu_list)) {
    seuobj <- seu_list[[i]]
    sample_id <- seuobj@meta.data$sample_id[1]
    if (verbose) message(paste("Start processing:", sample_id))

    # Normalisation and scaling
    seuobj <- NormalizeData(seuobj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    seuobj <- FindVariableFeatures(object = seuobj, nfeatures = 3000, verbose = FALSE, selection.method = "vst")
    # Optionally exclude some variable genes
    var_genes <- VariableFeatures(seuobj)
    var_genes <- var_genes[!var_genes %in% var_genes_to_exclude]
    VariableFeatures(seuobj) <- var_genes
    seuobj <- ScaleData(seuobj, vars.to.regress = "percent.mt", verbose = FALSE)
    seuobj <- RunPCA(seuobj, npcs = pca_dims, verbose = FALSE)
    seuobj <- RunUMAP(seuobj,
      reduction = "pca", dims = 1:pca_dims, n.neighbors = umap_neighbors, min.dist = umap_min_dist, n.epochs = umap_epochs,
      local.connectivity = umap_local_connectivity, seed.use = seed
    )
    seuobj <- FindNeighbors(seuobj, dims = 1:pca_dims, verbose = FALSE)
    seuobj <- FindClusters(seuobj, resolution = cluster_resolution, algorithm = cluster_algorithm, verbose = FALSE)

    # Store in the list
    seu_list[[i]] <- seuobj
  }

  return(seu_list)
}


basicplot_clusters <- function(seu_list, group.by = NULL, fig_dir = "Figures/basicplots/", output_name = "_plots.pdf") {
  # Create output directory if it doesn't exist
  if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

  for (i in seq_along(seu_list)) {
    seuobj <- seu_list[[i]]
    sample_id <- seuobj@meta.data$sample_id[1]
    sample_id <- sub("^(E[0-9]+_[A-Za-z]+).*", "\\1", sample_id)
    output_pdf <- file.path(fig_dir, paste0(sample_id, output_name))

    message("Generating plots for sample: ", sample_id)

    # Use active identity if group.by is NULL
    p1 <- DimPlot(
      seuobj,
      reduction = "umap",
      group.by = group.by,
      pt.size = 0.5,
      label = TRUE,
      label.size = 4
    ) + ggtitle(paste("UMAP -", sample_id)) +
      theme(plot.title = element_text(hjust = 0.5))

    p2 <- VlnPlot(
      seuobj,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
      pt.size = 0.2,
      group.by = group.by
    ) + ggtitle(paste("QC Metrics -", sample_id)) +
      theme(plot.title = element_text(hjust = 0.5))

    # Save plots to PDF
    pdf(output_pdf, width = 12, height = 10)
    print(plot_grid(p1, p2, ncol = 1))
    dev.off()

    message("Plots saved to: ", output_pdf)
  }
}

save_all_clusters_to_one_pdf <- function(
  seurat_list,
  out_pdf = "All_Datasets_Cluster_Highlight_Plots.pdf",
  highlight_color = "gold",
  background_color = "lightgray",
  reduction = "umap",
  pt.size = 1
) {
  # Start multi-page PDF
  pdf(out_pdf, width = 6, height = 5)

  # Loop through datasets and clusters
  for (dataset_name in names(seurat_list)) {
    seurat_obj <- seurat_list[[dataset_name]]
    cluster_ids <- sort(unique(Idents(seurat_obj)))

    for (cluster in cluster_ids) {
      plot_obj <- Cluster_Highlight_Plot(
        seurat_object = seurat_obj,
        cluster_name = as.character(cluster),
        highlight_color = highlight_color,
        background_color = background_color,
        reduction = reduction,
        pt.size = pt.size
      ) + ggplot2::ggtitle(paste0(dataset_name, "Cluster ", cluster))

      print(plot_obj)
    }
  }

  dev.off()
  message("Combined PDF saved")
}

batch_process_seurat_objects_SCTransform <- function(
  seu_list,
  percent_mt_regress = "percent.mt",
  var_genes_to_exclude = c("Xist", "Tsix", "Ddx3y", "Uty", "Eif2s3y", "Kdm5d", "Gm29650", "Malat1", "Gm42418", "AY036118", "Gm13305"),
  pca_dims = 50,
  umap_neighbors = 55L,
  umap_min_dist = 0.40,
  umap_epochs = 1000,
  umap_local_connectivity = 1L,
  seed = 11,
  cluster_resolution = 1.8,
  cluster_algorithm = 4,
  verbose = TRUE
) { # Initialize an empty list to store Seurat objects

  for (i in seq_along(seu_list)) {
    seuobj <- seu_list[[i]]
    sample_id <- seuobj@meta.data$sample_id[1]
    sample_id <- sub("^(E[0-9]+_[A-Za-z]+).*", "\\1", sample_id)
    if (verbose) message(paste("Start processing:", sample_id))

    # SCTransform
    message(paste0("Performing SCTransform for ", sample_id))
    seuobj <- SCTransform(seuobj, vars.to.regress = percent_mt_regress, verbose = verbose)

    # Identify variable genes and exclude specified ones
    message("Identifying and filtering variable genes...")
    var_genes <- VariableFeatures(object = seuobj)
    var_genes <- var_genes[!var_genes %in% var_genes_to_exclude]
    # VariableFeatures(seuobj) <- var_genes

    # PCA
    message("Running PCA...")
    seuobj <- RunPCA(seuobj, verbose = verbose, seed.use = seed, features = var_genes)

    # UMAP
    message("Running UMAP...")
    seuobj <- RunUMAP(seuobj,
      reduction = "pca", dims = 1:pca_dims, n.neighbors = umap_neighbors, min.dist = umap_min_dist,
      n.epochs = umap_epochs, local.connectivity = umap_local_connectivity, seed.use = seed
    )

    # Find neighbors
    message("Finding neighbors...")
    seuobj <- FindNeighbors(seuobj, dims = 1:pca_dims, verbose = FALSE)

    # Find clusters
    message("Finding clusters...")
    py_module_available(module = "leidenalg")
    seuobj <- FindClusters(seuobj, verbose = FALSE, resolution = cluster_resolution, algorithm = cluster_algorithm)

    # Plot UMAP
    message("Plotting UMAP...")
    DimPlot(seuobj, label = TRUE, reduction = "umap")

    # Store in Seurat objects in a list
    seu_list[[sample_id]] <- seuobj
  }

  print("All object processed successfully!")
  return(seu_list)
}

# Function to create EnhancedVolcano plot
createVolcanoPlot <- function(data, title, interest_genes) {
  EnhancedVolcano(data,
    lab = data$gene,
    x = "avg_logFC",
    y = "p_val_adj",
    FCcutoff = 0.25,
    pCutoff = 0.05,
    title = title,
    ylim = c(0, 28),
    legendPosition = "none",
    drawConnectors = TRUE,
    boxedLabels = TRUE,
    legendLabels = c("NS.", "LogFC", "p-value", "p-value and LogFC"),
    selectLab = interest_genes,
    max.overlaps = Inf,
    border = "full",
    labSize = 3,
    xlab = "log2(Fold Change)",
    ylab = expression("FDR adjusted -log"["10"] * italic("P"))
  )
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
