#!/usr/bin/env Rscript
# Step 2: Differential expression analysis of single-cell proteomics
#
# Differential protein analysis using DEP package.
# Treat each single cell as a replicate sample.
#
# Usage:
#   Rscript step2_de_analysis.R --config config.yaml
#   Rscript step2_de_analysis.R --config config.yaml --input preprocessed.tsv --exp_design exp_design.csv --output_dir results/step2

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})

# ============================================
# Parameter parsing
# ============================================
option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = "config.yaml",
              help = "Configuration file path (YAML)"),
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Preprocessed TSV file path (override config)"),
  make_option(c("-e", "--exp_design"), type = "character", default = NULL,
              help = "Experimental design CSV file path (override config)"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "Output directory (override config))")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ============================================
# Load configuration
# ============================================
cfg <- yaml.load_file(opt$config)

input_file <- if (!is.null(opt$input)) opt$input else {
  if (!is.null(cfg$de_analysis$input_file)) {
    cfg$de_analysis$input_file
  } else {
    file.path(cfg$output$dir, paste0(cfg$output$prefix, "_preprocessed.tsv"))
  }
}
exp_design_file <- if (!is.null(opt$exp_design)) opt$exp_design else cfg$input$exp_design
output_dir <- if (!is.null(opt$output_dir)) opt$output_dir else file.path(cfg$output$dir, "step2")

de_cfg <- cfg$de_analysis

# ============================================
# Check input files
# ============================================
if (!file.exists(input_file)) {
  stop(paste0("Error: Preprocessed file does not exist: ", input_file))
}
if (!file.exists(exp_design_file)) {
  stop(paste0("Error: Experimental design file does not exist: ", exp_design_file, "\nPlease run first generate_exp_design.py"))
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat(paste0(rep("=", 60), collapse = ""), "\n")
cat("Step 2: Differential expression analysis (DEP)\n")
cat(paste0(rep("=", 60), collapse = ""), "\n")
cat("Input data:", input_file, "\n")
cat("Experimental design:", exp_design_file, "\n")
cat("Output directory:", output_dir, "\n\n")

# ============================================
# Load R packages
# ============================================
suppressPackageStartupMessages({
  library(DEP)
  library(dplyr)
  library(RColorBrewer)
  library(SummarizedExperiment)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(umap)
  library(tibble)
})

# ============================================
# Reading data
# ============================================
cat("[1] Reading data...\n")
data <- fread(input_file, sep = "\t", header = TRUE, quote = "")
experiment_design <- read.csv(exp_design_file, sep = ",", header = TRUE)

# Get sample column indices (exclude PG.ProteinGroups and PG.Genes)
data_colums <- grep("_", colnames(data))

cat("  Data dimensions:", nrow(data), "rows x", ncol(data), "columns\n")
cat("  Number of samples:", length(data_colums), "\n")

# Convert NaN to NA (VSN does not accept NaN, but accepts NA)
for (col in data_colums) {
  data[[col]][is.nan(data[[col]])] <- NA
}
nan_count <- sum(is.nan(as.matrix(data[, ..data_colums])))
cat("  NaN→NA Conversion complete (remaining NaN:", nan_count, ")\n")

# ============================================
# Protein unique identifier
# ============================================
cat("\n[2] GenerateProtein unique identifier...\n")
data_unique <- make_unique(data, "PG.ProteinGroups", "PG.Genes", delim = ";")
data_se <- make_se(data_unique, data_colums, experiment_design)

cat("  SummarizedExperiment dimensions:", dim(data_se), "\n")

# ============================================
# Filter missing values
# ============================================
cat("\n[3] Filter missing values...\n")
thr <- de_cfg$filter_missval
if (!is.null(thr)) {
  data_se <- filter_missval(data_se, thr = thr)
  cat("  filterafterdimensions:", dim(data_se), "\n")
} else {
  cat("  Skip filtering (filter_missval = null)\n")
  cat("  Data dimensions:", dim(data_se), "\n")
}

# ============================================
# Data normalization
# ============================================
cat("\n[4] VSN Normalization...\n")
data_norm <- normalize_vsn(data_se)

# After normalize_vsn, assay row names are lost, need to recover using name column from rowData
# Must simultaneously set SE row names and assay row names, otherwise limma differential analysis results all Inf
if (is.null(rownames(assay(data_norm))) || all(rownames(assay(data_norm)) == "")) {
  protein_names <- make.names(rowData(data_norm)$name, unique = TRUE)
  rownames(data_norm) <- protein_names
  rownames(assay(data_norm)) <- protein_names
  cat("  Recovered assay rowsname:", length(protein_names), "proteins\n")
}

pdf(file.path(output_dir, "2norm_plot.pdf"), width = 6, height = max(6, ncol(data_norm) * 0.3))
plot_normalization(data_se, data_norm)
dev.off()
cat("  Saved: 2norm_plot.pdf\n")

# ============================================
# Missing value visualization and imputation
# ============================================
cat("\n[5] Missing value visualization and imputation...\n")
pdf(file.path(output_dir, "3missval_plot.pdf"), width = 12, height = 12)
plot_missval(data_norm)
plot_detect(data_norm)
dev.off()

impute_method <- de_cfg$imputation_method

if (impute_method == "by_group") {
  # ========================================
  # Impute by group minimum (recommended for single cell)
  # ========================================
  cat("  Using group minimum imputation method...\n")
  mat <- assay(data_norm)
  conditions <- colData(data_norm)$condition
  unique_groups <- unique(conditions)

  mat_imputed <- mat
  for (i in 1:nrow(mat)) {
    for (grp in unique_groups) {
      group_cols <- which(conditions == grp)
      current_protein_group_data <- mat[i, group_cols]

      if (all(is.na(current_protein_group_data))) {
        # This protein in this group is all NA -> fill 0
        mat_imputed[i, group_cols] <- 0
      } else if (any(is.na(current_protein_group_data))) {
        # This protein has partial missing in this group -> take group minimum value
        specific_min <- min(current_protein_group_data, na.rm = TRUE)
        na_indices <- is.na(current_protein_group_data)
        current_protein_group_data[na_indices] <- specific_min
        mat_imputed[i, group_cols] <- current_protein_group_data
      }
    }
  }
  data_imp <- data_norm
  assay(data_imp) <- mat_imputed

} else {
  # ========================================
  # Manually implement imputation method (avoid DEP::impute incompatibility with new SummarizedExperiment)
  # ========================================
  mat <- assay(data_norm)

  if (impute_method == "zero") {
    cat("  Using zero imputation method...\n")
    mat[is.na(mat)] <- 0
  } else if (impute_method == "MinDet") {
    cat("  Use MinDet imputemethod...\n")
    # DEP original MinDet (imputeLCMD::impute.MinDet):
    # Per sample (column), impute NA with q-th quantile of each sample (default 1%)
    # Suitable for left-censored missing values (below detection limit)
    q <- de_cfg$min_det_q %||% 0.01
    low_quantile <- apply(mat, 2, quantile, prob = q, na.rm = TRUE)
    for (i in 1:ncol(mat)) {
      na_idx <- is.na(mat[, i])
      if (any(na_idx)) {
        mat[na_idx, i] <- low_quantile[i]
      }
    }
    cat("  Use q =", q, "quantileimpute\n")
  } else if (impute_method == "MinProb") {
    cat("  Use MinProb imputemethod...\n")
    # DEP original MinProb (imputeLCMD::impute.MinProb):
    # Per sample (column), draw from normal distribution N(q-quantile - tune*sd, sd)
    q <- de_cfg$min_prob_q %||% 0.01
    tune <- de_cfg$min_prob_tune %||% 1
    low_quantile <- apply(mat, 2, quantile, prob = q, na.rm = TRUE)
    col_sds <- apply(mat, 2, sd, na.rm = TRUE)
    col_sds[is.na(col_sds) | col_sds == 0] <- 1
    for (i in 1:ncol(mat)) {
      na_idx <- is.na(mat[, i])
      if (any(na_idx)) {
        mat[na_idx, i] <- stats::rnorm(sum(na_idx),
                                       mean = low_quantile[i] - tune * col_sds[i],
                                       sd = col_sds[i])
      }
    }
  } else if (impute_method == "knn") {
    cat("  Use KNN imputemethod...\n")
    if (!requireNamespace("impute", quietly = TRUE)) {
      stop("KNN imputation needs 'impute' package, please install: BiocManager::install('impute')")
    }
    mat <- impute::impute.knn(as.matrix(mat))$data
  } else {
    stop(paste0("notsupportofimputemethod: ", impute_method))
  }

  data_imp <- data_norm
  assay(data_imp, withDimnames = FALSE) <- mat
  # After impute, reset assay row names (ensure limma differential analysis works)
  rownames(assay(data_imp)) <- rownames(data_norm)
}

# ============================================
# Differential analysis
# ============================================
cat("\n[6] Differential analysis...\n")
control_name <- de_cfg$control_name

data_diff <- test_diff(data_imp, type = "control", control = control_name)

alpha <- de_cfg$alpha
lfc_thr <- de_cfg$lfc_threshold

# Manually implement add_rejections logic (add_rejections has bug in new SummarizedExperiment)：
# Criteria: p.adj <= alpha and |diff| >= log2(lfc)(with DEP::add_rejections consistent)
dep <- data_diff
rd <- rowData(dep)
padj_cols <- grep("_p\\.adj$", colnames(rd), value = TRUE)
for (c_col in padj_cols) {
  base <- sub("_p\\.adj$", "", c_col)
  diff_col <- paste0(base, "_diff")
  sig_col <- paste0(base, "_significant")
  if (diff_col %in% colnames(rd)) {
    rd[[sig_col]] <- (rd[[c_col]] <= alpha) & (abs(rd[[diff_col]]) >= log2(lfc_thr))
  }
}
rd$significant <- rowSums(as.matrix(rd[, paste0(sub("_p\\.adj$", "", padj_cols), "_significant"), drop = FALSE]), na.rm = TRUE) > 0
rowData(dep) <- rd

n_sig <- sum(rd$significant, na.rm = TRUE)
cat("  significantdifferentialproteincount:", n_sig, "\n")

# ============================================
# UMAP visualization
# ============================================
cat("\n[7] UMAP visualization...\n")
expr_matrix <- assay(data_imp)
set.seed(123)

n_neighbors <- de_cfg$umap_n_neighbors %||% 5
samples_as_rows <- t(expr_matrix)
umap_result <- umap(samples_as_rows, n_components = 2, n_neighbors = n_neighbors)

group_ids <- colnames(expr_matrix)
condition <- colData(data_imp)$condition

umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = condition,
  row.names = group_ids,
  stringsAsFactors = FALSE
)

pdf(file.path(output_dir, "4umap_plot.pdf"), width = 8, height = 6)
p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "UMAP of Protein Expression", x = "UMAP1", y = "UMAP2") +
  theme_minimal()
print(p)

# UMAP with labels
umap_df <- rownames_to_column(umap_df, var = "SampleID")
p2 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_text(aes(label = SampleID), size = 2, vjust = -0.5, hjust = 0.5,
            color = "black", alpha = 0.6) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "UMAP (with labels)", x = "UMAP1", y = "UMAP2") +
  theme_minimal()
print(p2)
dev.off()
cat("  Saved: 4umap_plot.pdf\n")

# ============================================
# PCA, correlation, heatmap, volcano plot
# ============================================
cat("\n[8] PCA、heatmap、volcano plot...\n")

n_features <- de_cfg$pca_n_features
if (is.null(n_features)) n_features <- nrow(data_imp)

pdf(file.path(output_dir, "5pca_plot.pdf"), width = 8, height = 6)
pcaplot <- plot_pca(data_imp, x = 1, y = 2, n = n_features, point_size = 3, label = FALSE, indicate = "condition")
print(pcaplot)
dev.off()
cat("  Saved: 5pca_plot.pdf\n")

# Protein identification count
tryCatch({
  pdf(file.path(output_dir, "1se_plot.pdf"), width = 18, height = 6)
  p1 <- plot_frequency(data_se)
  print(p1)

  p2 <- DEP2::plot_numbers(data_se)
  print(p2)

  dev.off()
  cat("  Saved: 1se_plot.pdf\n")
}, error = function(e) {
  cat("  Warning: SE plot failed:", e$message, "\n")
  if (exists("pdf")) dev.off()
})

pdf(file.path(output_dir, "6heatmap_volcano.pdf"), width = 12, height = 12)

# Attempt to plot heatmap and volcano plot
tryCatch({
  # Correlation plot
  plot_cor(dep, significant = FALSE, lower = -1, upper = 1, pal = "GnBu")

  # heatmap
  plot_heatmap(dep, type = "centered", kmeans = TRUE, k = 6, col_limit = 4,
               show_row_names = FALSE, indicate = "condition")
  kmeanscluster <- plot_heatmap(dep, type = "centered", kmeans = TRUE, k = 6, col_limit = 4,
                                 show_row_names = FALSE, indicate = "condition", plot = FALSE)
  write.csv(kmeanscluster, file.path(output_dir, "6kmeanscluster.csv"))

  # volcano plot
  contrasts <- colnames(assay(dep))[grep("contrast", colnames(assay(dep)))]
  for (con in contrasts) {
    cat("  volcano plot contrast:", con, "\n")
    tryCatch({
      DEP2::plot_volcano(dep, contrast = con, label_size = 2, add_names = TRUE)
    }, error = function(e) {
      cat("    Warning:", e$message, "\n")
    })
  }
}, error = function(e) {
  cat("  Warning: heatmap/volcano plot generation failed:", e$message, "\n")
})

dev.off()
cat("  Saved: 6heatmap_volcano.pdf\n")

# ============================================
# Extract results
# ============================================
cat("\n[9] Extract results...\n")
data_results <- get_results(dep)
write.csv(data_results, file.path(output_dir, "7data_results.csv"), row.names = FALSE)

# Save R objects
save(data_se, data_norm, data_imp, data_diff, dep,
     file = file.path(output_dir, "data.RData"))

cat("\n", paste0(rep("=", 60), collapse = ""), "\n")
cat("Step 2 complete!\n")
cat("resultsSavein:", output_dir, "\n")

# List all generated files
cat("\nOutput files:\n")
output_files <- list.files(output_dir, pattern = "\\.(pdf|csv|tsv|txt)$", full.names = FALSE)
for (f in output_files) {
  file_path <- file.path(output_dir, f)
  file_size <- file.info(file_path)$size
  if (file_size < 1024) {
    size_str <- paste0(file_size, " B")
  } else if (file_size < 1024^2) {
    size_str <- paste0(round(file_size / 1024, 1), " KB")
  } else {
    size_str <- paste0(round(file_size / 1024^2, 1), " MB")
  }
  cat("  -", f, "(", size_str, ")\n")
}
cat(paste0(rep("=", 60), collapse = ""), "\n")
