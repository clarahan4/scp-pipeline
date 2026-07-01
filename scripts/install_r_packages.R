#!/usr/bin/env Rscript
# Install required R packages for this pipeline
#
# Usage:
#   Rscript install_r_packages.R
#
# - Bioconductor packages via BiocManager
# - CRAN packages via install.packages
# - DEP2 from GitHub

# BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

cat("============================================\n")
cat("  Installing R packages for SCP pipeline\n")
cat("============================================\n\n")

# Bioconductor packages
bioc_pkgs <- c(
  "DEP",
  "SummarizedExperiment",
  "Seurat",
  "monocle3"
)

cat("[1] Installing Bioconductor packages...\n")
cat("    ", paste(bioc_pkgs, collapse = ", "), "\n")
BiocManager::install(bioc_pkgs, update = FALSE, ask = FALSE, checkBuilt = FALSE)

# CRAN packages
cran_pkgs <- c(
  "yaml", "optparse", "dplyr", "ggplot2", "patchwork",
  "RColorBrewer", "rmarkdown", "pheatmap", "data.table",
  "tibble", "umap", "limma", "Matrix", "tidyverse", "ggsci"
)

cat("\n[2] Installing CRAN packages...\n")
cat("    ", paste(cran_pkgs, collapse = ", "), "\n")
install.packages(cran_pkgs, repos = "https://cloud.r-project.org", checkBuilt = FALSE)

# DEP2 from GitHub
cat("\n[3] Installing DEP2 (GitHub)...\n")
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
}
devtools::install_github("cilikes/DEP2", build_vignettes = FALSE, force = FALSE)

cat("\n============================================\n")
cat("  Done!\n")
cat("============================================\n")

# Verify
cat("\n[4] Verifying...\n")
all_pkgs <- c(bioc_pkgs, cran_pkgs, "DEP2", "devtools")
for (pkg in all_pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    version <- tryCatch(
      as.character(packageVersion(pkg)),
      error = function(e) "unknown"
    )
    cat(sprintf("  OK  %-30s %s\n", pkg, version))
  } else {
    cat(sprintf("  X   %-30s not installed\n", pkg))
  }
}
