#!/usr/bin/env Rscript
# Step 3: Single-cell Proteomics Clustering Analysis (Seurat)
#
# Perform single-cell proteomics dimensionality reduction, clustering, visualization, and trajectory analysis using Seurat.
#
# Usage:
#   Rscript step3_cluster_seurat.R --config config.yaml
#   Rscript step3_cluster_seurat.R --config config.yaml --input preprocessed.tsv --output_dir results/step3_seurat

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})

# ============================================
# Argument Parsing
# ============================================
option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = "config.yaml",
              help = "Configuration file path (YAML)"),
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to preprocessed TSV file (overrides config)"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "Output directory (overrides config)")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ============================================
# Load Configuration
# ============================================
cfg <- yaml.load_file(opt$config)

input_file <- if (!is.null(opt$input)) opt$input else {
  if (!is.null(cfg$seurat$input_file)) {
    cfg$seurat$input_file
  } else {
    file.path(cfg$output$dir, paste0(cfg$output$prefix, "_preprocessed.tsv"))
  }
}
output_dir <- if (!is.null(opt$output_dir)) opt$output_dir else file.path(cfg$output$dir, "step3_seurat")
seurat_cfg <- cfg$seurat

if (!file.exists(input_file)) {
  stop(paste0("Error: Preprocessed file does not exist: ", input_file))
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat(paste0(rep("=", 60), collapse = ""), "\n")
cat("Step 3: Seurat Single-cell Clustering Analysis\n")
cat(paste0(rep("=", 60), collapse = ""), "\n")
cat("Input data:", input_file, "\n")
cat("Output directory:", output_dir, "\n\n")

# ============================================
# Load R Packages
# ============================================
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(RColorBrewer)
  library(limma)
})

# ============================================
# Read Data
# ============================================
cat("[1] Reading data...\n")
raw <- read.table(input_file, sep = "\t", header = TRUE, quote = "", check.names = FALSE)

# Create row names: ProteinID_GeneID
row_names <- paste(raw$PG.ProteinGroups, raw$PG.Genes, sep = "_")

# Extract expression matrix (remove first two metadata columns)
expr <- raw %>% select(-PG.ProteinGroups, -PG.Genes) %>% as.matrix()
rownames(expr) <- row_names

cat("  Expression matrix dimensions:", nrow(expr), "proteins x", ncol(expr), "cells\n")

# ============================================
# Missing Value Imputation
# ============================================
cat("\n[2] Imputing missing values...\n")
na_method <- seurat_cfg$na_fill_method %||% "zero"

expr_original <- expr  # Keep original matrix for QC

if (na_method == "zero") {
  expr[is.na(expr) | is.nan(expr)] <- 0
  cat("  Using 0 to fill missing values\n")
} else if (na_method == "min") {
  expr <- apply(expr, 1, function(row) {
    row_min <- min(row, na.rm = TRUE)
    row[is.na(row) | is.nan(row)] <- row_min
    return(row)
  }) %>% t()
  cat("  Using row minimum to fill missing values\n")
}

# ============================================
# Create Seurat Object
# ============================================
cat("\n[3] Creating Seurat object...\n")
cell_names <- colnames(expr)

# Parse condition from column names
rename_pattern <- cfg$preprocess$rename_pattern
if (!is.null(rename_pattern) && nchar(rename_pattern) > 0) {
  # Try to parse from simplified column names (e.g., control_1 → control)
  group_info <- sapply(cell_names, function(name) {
    parts <- strsplit(name, "_")[[1]]
    return(parts[1])
  })
} else {
  group_info <- rep("unknown", length(cell_names))
}

min_cells <- seurat_cfg$min_cells %||% 1
min_features <- seurat_cfg$min_features %||% 1

obj <- CreateSeuratObject(counts = expr, project = "SCP",
                         min.cells = min_cells, min.features = min_features)
obj$group <- group_info
obj$cell_type <- "SCP"

cat("  Seurat object:", ncol(obj), "cells x", nrow(obj), "proteins\n")
cat("  Group statistics:\n")
print(table(obj$group))

# ============================================
# QC: Protein Count
# ============================================
cat("\n[4] Protein count QC...\n")
protein_counts <- apply(expr_original, 2, function(x) {
  sum(!is.na(x) & !is.nan(x))
})

obj$protein_count <- protein_counts[colnames(obj)]

plot_data <- data.frame(
  Cell = names(protein_counts[colnames(obj)]),
  Protein_Count = as.numeric(protein_counts[colnames(obj)]),
  Group = group_info,
  stringsAsFactors = FALSE
)

pdf(file.path(output_dir, "1protein_count.pdf"), width = max(8, ncol(obj) * 0.15), height = 5)
ggplot(plot_data, aes(x = Cell, y = Protein_Count, fill = Group)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Protein Count per Cell", x = "Cell", y = "Number of Proteins") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5))
dev.off()
cat("  Saved: 1protein_count.pdf\n")

# VlnPlot: Seurat 5.2's VlnPlot internally still uses the defunct slot parameter,
# here we use ggplot2 to manually draw violin plots (based on nCount_RNA = total protein amount)
nCount_RNA <- colSums(expr_original, na.rm = TRUE)[colnames(obj)]
vln_data <- data.frame(
  nCount_RNA = as.numeric(nCount_RNA),
  group = obj$group,
  stringsAsFactors = FALSE
)
pdf(file.path(output_dir, "1vlnplot.pdf"), width = 8, height = 4)
print(ggplot(vln_data, aes(x = group, y = nCount_RNA, fill = group)) +
        geom_violin(trim = FALSE, alpha = 0.7) +
        geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
        scale_fill_brewer(palette = "Set1") +
        labs(title = "nCount_RNA by Group", x = "Group", y = "nCount_RNA") +
        theme_minimal())
dev.off()
cat("  Saved: 1vlnplot.pdf\n")

cat("  Protein count statistics summary:\n")
print(summary(protein_counts))

# ============================================
# Standard Analysis Pipeline
# ============================================
cat("\n[5] Standard analysis pipeline...\n")
obj_ana <- obj

norm_method <- seurat_cfg$normalization %||% "LogNormalize"
scale_factor <- seurat_cfg$scale_factor %||% 1000
n_var_features <- seurat_cfg$n_variable_features %||% 800
pca_dims <- seurat_cfg$pca_dims %||% 10
n_neighbors <- seurat_cfg$n_neighbors %||% 20
min_dist <- seurat_cfg$min_dist %||% 0.05
resolution <- seurat_cfg$resolution %||% 1.0

if (norm_method == "SCTransform") {
  cat("  Using SCTransform...\n")
  obj_ana <- SCTransform(obj_ana, variable.features.n = n_var_features, verbose = FALSE)
} else {
  cat("  Using LogNormalize (scale_factor =", scale_factor, ")...\n")
  obj_ana <- NormalizeData(obj_ana, normalization.method = "LogNormalize", scale.factor = scale_factor)
  # Merge layers (required for Seurat 5, to avoid slot defunct bug in subsequent VlnPlot/FindAllMarkers)
  obj_ana[["RNA"]] <- JoinLayers(obj_ana[["RNA"]])
  obj_ana <- FindVariableFeatures(obj_ana, selection.method = "vst", nfeatures = n_var_features)
  obj_ana <- ScaleData(obj_ana)
}

obj_ana <- RunPCA(obj_ana, verbose = FALSE)

# ElbowPlot
pdf(file.path(output_dir, "2elbow_plot.pdf"), width = 6, height = 4)
print(ElbowPlot(obj_ana, ndims = 30))
dev.off()
cat("  Saved: 2elbow_plot.pdf\n")

cat("  UMAP (dims = 1:", pca_dims, ", n_neighbors =", n_neighbors, ", min_dist =", min_dist, ")...\n")
obj_ana <- RunUMAP(obj_ana, dims = 1:pca_dims, reduction = "pca",
                   n.neighbors = n_neighbors, min.dist = min_dist)
obj_ana <- FindNeighbors(obj_ana, dims = 1:pca_dims, reduction = "pca")
obj_ana <- FindClusters(obj_ana, resolution = resolution)

n_clusters <- length(unique(obj_ana$seurat_clusters))
cat("  Identified", n_clusters, "clusters\n")

# ============================================
# UMAP Visualization
# ============================================
cat("\n[6] UMAP visualization...\n")

# Automatic color assignment
n_cl <- n_clusters
if (n_cl <= 4) {
  cluster_colors <- c("#F8766D", "#7CAE01", "#01BFC4", "#C77CFF")
} else {
  cluster_colors <- colorRampPalette(brewer.pal(8, "Set1"))(n_cl)
}

pdf(file.path(output_dir, "3umap.pdf"), width = 8, height = 5)
print(DimPlot(obj_ana, reduction = "umap", group.by = "group", pt.size = 0.5) +
        ggtitle("Colored by Group"))
print(DimPlot(obj_ana, reduction = "umap", group.by = "seurat_clusters",
              pt.size = 0.5, cols = cluster_colors) +
        ggtitle("Clusters"))
dev.off()
cat("  Saved: 3umap.pdf\n")

# ============================================
# Marker Protein Analysis
# ============================================
cat("\n[7] Marker protein analysis...\n")
marker_test <- seurat_cfg$marker_test %||% "wilcox"
marker_logfc <- seurat_cfg$marker_logfc_threshold %||% 0.25
marker_min_pct <- seurat_cfg$marker_min_pct %||% 0.25
top_n <- seurat_cfg$top_n_markers %||% 10

obj_markers <- FindAllMarkers(obj_ana,
                              logfc.threshold = marker_logfc,
                              test.use = marker_test,
                              min.pct = marker_min_pct,
                              assay = "RNA",
                              only.pos = TRUE,
                              verbose = TRUE,
                              layer = "data")

# Error tolerance: FindAllMarkers may return empty in single-cell proteomics (No DE genes)
if (is.null(obj_markers) || nrow(obj_markers) == 0) {
  cat("  Warning: No marker proteins detected (criteria too strict), relaxing thresholds and retrying...\n")
  obj_markers <- FindAllMarkers(obj_ana,
                                logfc.threshold = 0,
                                test.use = marker_test,
                                min.pct = 0,
                                assay = "RNA",
                                only.pos = TRUE,
                                verbose = FALSE,
                                layer = "data")
}

write.table(obj_markers, file = file.path(output_dir, "4markers.csv"),
            quote = FALSE, col.name = TRUE, row.names = TRUE, sep = "\t")
cat("  Saved: 4markers.csv, total", nrow(obj_markers), "markers\n")

# Top N markers
if (nrow(obj_markers) > 0) {
  top_markers <- obj_markers %>% group_by(cluster) %>% top_n(n = top_n, wt = avg_log2FC)
  write.csv(top_markers, file = file.path(output_dir, "5top_markers.csv"), row.names = FALSE)
} else {
  top_markers <- data.frame()
  cat("  Warning: No markers, skipping top markers and dotplot\n")
}

# DotPlot and Heatmap (skip when markers are empty to avoid errors)
if (nrow(top_markers) > 0) {
  tryCatch({
    pdf(file.path(output_dir, "6markers_dotplot_heatmap.pdf"), width = 12, height = 8)
    print(DotPlot(obj_ana, cols = c("lightgrey", "#E64B35FF"),
                  features = unique(top_markers$gene)))
    print(DoHeatmap(object = obj_ana, features = top_markers$gene, layer = "data"))
    dev.off()
    cat("  Saved: 6markers_dotplot_heatmap.pdf\n")
  }, error = function(e) {
    cat("  Warning: DotPlot/Heatmap plotting failed:", e$message, "\n")
    dev.off()
  })
} else {
  cat("  Skip DotPlot/Heatmap (no markers)\n")
}

# ============================================
# Save Results
# ============================================
cat("\n[9] Saving results...\n")
save(obj_ana, file = file.path(output_dir, "seurat_result.RData"))
write.csv(obj_ana@meta.data, file = file.path(output_dir, "metadata.csv"), row.names = TRUE)

cat("\n", paste0(rep("=", 60), collapse = ""), "\n")
cat("Step 3 completed!\n")
cat("Results saved in:", output_dir, "\n")
cat(paste0(rep("=", 60), collapse = ""), "\n")
