#!/usr/bin/env Rscript
# Old environment R package installation (install strictly in scp-pipeline-old environment, R 4.4.1 + Bioc 3.19)
# Used to replicate original DEP analysis (DEP 1.26 original, no assays<- compatibility issues)
#
# Usage (must use conda run to specify environment, avoid installing to base):
#   conda run -n scp-pipeline-old Rscript scripts/install_r_packages_old.R

# Confirm installing to environment library
cat("R version:", R.version.string, "\n")
cat(".libPaths():\n")
print(.libPaths())
env_lib <- .libPaths()[1]
cat("Installing to:", env_lib, "\n\n")

# BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = env_lib)
}

cat("============================================\n")
cat("  Installing old R packages (Bioc 3.19, DEP 1.26 original)\n")
cat("============================================\n\n")

# Bioconductor 3.19: DEP 1.26 (includes SummarizedExperiment 1.34 old version / limma dependencies)
# force=TRUE ensures compilation from source (avoid pulling R4.5 precompiled binaries)
cat("[1] Installing DEP (Bioc 3.19, compile from source)...\n")
BiocManager::install("DEP", version = "3.19", update = FALSE, ask = FALSE,
                     checkBuilt = FALSE, force = TRUE, lib = env_lib)

cat("\n============================================\n")
cat("  Verify versions\n")
cat("============================================\n")
verify <- c("DEP", "SummarizedExperiment", "limma", "imputeLCMD",
            "data.table", "ggplot2", "umap", "tibble")
for (pkg in verify) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    built <- utils::packageDescription(pkg)$Built
    cat(sprintf("  OK  %-25s %s  (Built: %s)\n",
                pkg, as.character(packageVersion(pkg)), substr(built, 1, 30)))
  } else {
    cat(sprintf("  X   %-25s not installed\n", pkg))
  }
}
