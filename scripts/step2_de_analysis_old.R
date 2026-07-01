#!/usr/bin/env Rscript
# Step 2 (old version environment): Differential expression analysis of single-cell proteomics data
#  scp-pipeline-old (R 4.3 + DEP 1.24 + imputeLCMD)
# Fully utilize the original DEP functions: impute() and add_rejections()
#
# Usage:
#   conda activate scp-pipeline-old
#   Rscript scripts/step2_de_analysis_old.R --config config.yaml

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})

# ============================================
# Parameter parsing
# ============================================
option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = "config.yaml",
              help = "YAML file"),
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "TSV file"),
  make_option(c("-e", "--exp_design"), type = "character", default = NULL,
              help = "CSV file"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "output dir")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ============================================
# Load configuration
# ============================================
cfg <- yaml.load_file(opt$config)

input_file <- if (!is.null(opt$input)) opt$input else {
  file.path(cfg$output$dir, paste0(cfg$output$prefix, "_preprocessed.tsv"))
}
exp_design_file <- if (!is.null(opt$exp_design)) opt$exp_design else cfg$input$exp_design
output_dir <- if (!is.null(opt$output_dir)) opt$output_dir else file.path(cfg$output$dir, "step2")

de_cfg <- cfg$de_analysis

if (!file.exists(input_file)) stop(paste0("Error: Preprocessed file does not exist: ", input_file))
if (!file.exists(exp_design_file)) stop(paste0("Error: Experimental design file does not exist:", exp_design_file))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat(paste0(rep("=", 60), collapse = ""), "\n")
cat("Step2: \n")
cat(paste0(rep("=", 60), collapse = ""), "\n")
cat("Input data:", input_file, "\n")
cat("Experiment file:", exp_design_file, "\n")
cat("Output dir:", output_dir, "\n")
cat("DEP version:", as.character(packageVersion("DEP")), "\n\n")

# ============================================
# Loading R packages
# ============================================
suppressPackageStartupMessages({
  library(DEP)
  library(SummarizedExperiment)
  library(data.table)
  library(ggplot2)
  library(umap)
  library(tibble)
})

# ============================================
# Reading data
# ============================================
cat("[1] Reading data...\n")
data <- fread(input_file, sep = "\t", header = TRUE, quote = "")
experiment_design <- read.csv(exp_design_file, sep = ",", header = TRUE)

data_colums <- grep("_", colnames(data))
cat("  data:", nrow(data), "row x", ncol(data), "column\n")
cat("  sample:", length(data_colums), "\n")

# NaN -> NA
for (col in data_colums) {
  data[[col]][is.nan(data[[col]])] <- NA
}
cat("  NaN->NA finished\n")

# ============================================
# Protein unique identifier
# ============================================
cat("\n[2] Protein unique identifier...\n")
data_unique <- make_unique(data, "PG.ProteinGroups", "PG.Genes", delim = ";")
data_se <- make_se(data_unique, data_colums, experiment_design)
cat("  SummarizedExperiment :", dim(data_se), "\n")

# ============================================
# Filter missing values
# ============================================
cat("\n[3] Filter missing values...\n")
thr <- de_cfg$filter_missval
if (!is.null(thr)) {
  data_se <- filter_missval(data_se, thr = thr)
  cat("  After filtering:", dim(data_se), "\n")
} else {
  cat("  (filter_missval = null)\n")
}

# ============================================
# VSN normalization
# ============================================
cat("\n[4] VSN Normalization...\n")
data_norm <- normalize_vsn(data_se)

# After normalize_vsn, assay row names are lost, need to recover using name column from rowData
# Must simultaneously set SE row names and assay row names, otherwise DEP::impute reports assays<- error
if (is.null(rownames(assay(data_norm))) || all(rownames(assay(data_norm)) == "")) {
  protein_names <- make.names(rowData(data_norm)$name, unique = TRUE)
  rownames(data_norm) <- protein_names
  rownames(assay(data_norm)) <- protein_names
  cat("  Recovered assay row names:", length(protein_names), "proteins\n")
}

pdf(file.path(output_dir, "2norm_plot.pdf"), width = 6, height = max(6, ncol(data_norm) * 0.3))
plot_normalization(data_se, data_norm)
dev.off()
cat("  Saved: 2norm_plot.pdf\n")

# ============================================
# Missing value visualization and imputation (DEP original impute)
# ============================================
cat("\n[5] Missing value visualization and imputation...\n")
pdf(file.path(output_dir, "3missval_plot.pdf"), width = 12, height = 12)
plot_missval(data_norm)
plot_detect(data_norm)
dev.off()

impute_method <- de_cfg$imputation_method

if (impute_method == "by_group") {
  # Impute by group minimum
  cat("  Using group minimum imputation method...\n")
  mat <- assay(data_norm)
  conditions <- colData(data_norm)$condition
  unique_groups <- unique(conditions)
  mat_imputed <- mat
  for (i in 1:nrow(mat)) {
    for (grp in unique_groups) {
      group_cols <- which(conditions == grp)
      current <- mat[i, group_cols]
      if (all(is.na(current))) {
        mat_imputed[i, group_cols] <- 0
      } else if (any(is.na(current))) {
        specific_min <- min(current, na.rm = TRUE)
        na_indices <- is.na(current)
        current[na_indices] <- specific_min
        mat_imputed[i, group_cols] <- current
      }
    }
  }
  data_imp <- data_norm
  assay(data_imp) <- mat_imputed
} else {
  # DEP original impute (MinDet/MinProb/knn/zero)
  cat("  Using DEP original impute, fun =", impute_method, "...\n")
  if (impute_method == "zero") {
    data_imp <- impute(data_norm, fun = "MinDet", q = 0.0001)
  } else {
    data_imp <- impute(data_norm, fun = impute_method)
  }
}

# ============================================
# Differential analysis (DEP original test_diff + add_rejections)
# ============================================
cat("\n[6] Differential analysis...\n")
control_name <- de_cfg$control_name

data_diff <- test_diff(data_imp, type = "control", control = control_name)

alpha <- de_cfg$alpha
lfc_thr <- de_cfg$lfc_threshold
dep <- add_rejections(data_diff, alpha = alpha, lfc = log2(lfc_thr))

n_sig <- rowData(dep)$significant %>% sum(na.rm = TRUE)
cat("  Significant differential proteins:", n_sig, "\n")

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
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = group_ids,
  stringsAsFactors = FALSE
)

pdf(file.path(output_dir, "4umap_plot.pdf"), width = 8, height = 6)
p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "UMAP of Protein Expression", x = "UMAP1", y = "UMAP2") +
  theme_minimal()
print(p)

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
# PCA, heatmap, volcano plot
# ============================================
cat("\n[8] PCA, heatmap, volcano plot...\n")

pdf(file.path(output_dir, "1se_plot.pdf"), width = 18, height = 6)
plot_frequency(data_se)
DEP2::plot_numbers(data_se)
dev.off()
cat("  Saved: 1se_plot.pdf\n")

n_features <- de_cfg$pca_n_features
if (is.null(n_features)) n_features <- nrow(data_imp)
pdf(file.path(output_dir, "5pca_plot.pdf"), width = 8, height = 6)
print(plot_pca(data_imp, x = 1, y = 2, n = n_features, point_size = 3, label = FALSE))
dev.off()
cat("  Saved: 5pca_plot.pdf\n")

pdf(file.path(output_dir, "6heatmap_volcano.pdf"), width = 12, height = 12)
tryCatch({
  plot_cor(dep, significant = FALSE, lower = -1, upper = 1, pal = "GnBu")
  plot_heatmap(dep, type = "centered", kmeans = TRUE, k = 6, col_limit = 4,
               show_row_names = FALSE, indicate = "condition")
  kmeanscluster <- plot_heatmap(dep, type = "centered", kmeans = TRUE, k = 6, col_limit = 4,
                                 show_row_names = FALSE, indicate = "condition", plot = FALSE)
  write.csv(kmeanscluster, file.path(output_dir, "6kmeanscluster.csv"))
  contrasts <- colnames(assay(dep))[grep("contrast", colnames(assay(dep)))]
  for (con in contrasts) {
    tryCatch(DEP2::plot_volcano(dep, contrast = con, label_size = 2, add_names = TRUE),
             error = function(e) cat("    Volcano plot warning:", e$message, "\n"))
  }
}, error = function(e) cat("  Warning: Heatmap/volcano plot failed:", e$message, "\n"))
dev.off()
cat("  Saved: 6heatmap_volcano.pdf\n")

# ============================================
# Extract results
# ============================================
cat("\n[9] Extract results...\n")
data_results <- get_results(dep)
write.csv(data_results, file.path(output_dir, "7data_results.csv"), row.names = FALSE)

save(data_se, data_norm, data_imp, data_diff, dep,
     file = file.path(output_dir, "data.RData"))

cat("\n", paste0(rep("=", 60), collapse = ""), "\n")
cat("Step 2 (old version environment) complete!\n")
cat("Significant differential proteins:", n_sig, "\n")
cat("Results saved at:", output_dir, "\n")
cat(paste0(rep("=", 60), collapse = ""), "\n")
