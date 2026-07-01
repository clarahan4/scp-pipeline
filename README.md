# Single-Cell Proteomics Analysis Pipeline

A config-driven pipeline for analyzing single-cell proteomics data from Spectronaut DIA output.

## Overview

```
Spectronaut raw data (TSV)
       │
       ▼
[Step 1] Data preprocessing ── Clean format, filter proteins, standardize sample names
       │
       ├──▶ [Step 2] Differential expression analysis (DEP/R) ── Identify differential proteins
       │
       ├──▶ [Step 3] Seurat clustering (R) ── Dimensionality reduction, clustering, UMAP
       │
       └──▶ [Step 4] Scanpy clustering (Python) ── Dimensionality reduction, clustering, Leiden, UMAP
              (Choose either Step 3 or Step 4)
```

## Quick Start

### Prerequisites

- Python >= 3.10
- R >= 4.3
- Conda or Mamba

### Installation

```bash
# 1. Clone repository
git clone https://github.com/clarahan4/scp-pipeline.git
cd scp-pipeline

# 2. Create Conda environment (recommended)
conda env create -f env/environment.yml
conda activate scp-pipeline

# 3. Install R packages
Rscript scripts/install_r_packages.R

# For older R version compatibility (R 4.4.1), use:
# conda env create -f env/environment_old.yml
# conda activate scp-pipeline-old
# Rscript scripts/install_r_packages_old.R
```

### Run the Pipeline

#### Method 1: Run complete workflow in one command (recommended)

```bash
# Run complete workflow (step1 → step2 → step3)
bash run_pipeline.sh config.yaml

# Use default config file
bash run_pipeline.sh
```

Which steps to run are determined by the `steps` configuration in `config.yaml`:

```yaml
steps:
  step1: true    # Preprocessing
  step2: true    # Differential expression analysis
  step3: true    # Seurat clustering
  step4: true   # Scanpy clustering
```

#### Method 2: Select which steps to run

Flexibly combine steps via command-line arguments (overrides config.yaml):

```bash
# Only run step1 + step2 (preprocessing + differential analysis)
bash run_pipeline.sh config.yaml --step1 --step2

# Only run step1 + step3 (preprocessing + Seurat clustering)
bash run_pipeline.sh config.yaml --step1 --step3

# Only run step1 + step4 (preprocessing + Scanpy clustering)
bash run_pipeline.sh config.yaml --step1 --step4

# Skip step1, only run step2 + step3 (use existing preprocessed data)
bash run_pipeline.sh config.yaml --step2 --step3

# Run all steps
bash run_pipeline.sh config.yaml --step1 --step2 --step3 --step4

# Only run clustering steps (step3 + step4)
bash run_pipeline.sh config.yaml --step3 --step4
```

### Using Custom Preprocessed Files

If you have your own preprocessed file and want to run Step 2/3/4 directly:

```bash
# Create a symbolic link to your custom file
mkdir -p results
ln -sf /path/to/your/preprocessed_file.tsv results/scp_analysis_preprocessed.tsv

# Then run any downstream steps
bash run_pipeline.sh config.yaml --step2 --step3
bash run_pipeline.sh config.yaml --step4
```

Note: Your custom file must follow the same format as Step 1 output (TSV with protein rows and sample columns, including `PG.ProteinGroups` and `PG.Genes` columns).

### Running Individual Steps

You can run any step individually (as long as its input file exists):

#### Run Step 2 individually (Differential expression analysis)

```bash
# Use configuration in config.yaml
Rscript scripts/step2_de_analysis.R --config config.yaml

# Manually specify input file (override config)
Rscript scripts/step2_de_analysis.R \
  --config config.yaml \
  --input my_custom_data.tsv \
  --output_dir results/step2
```

#### Run Step 3 individually (Seurat clustering)

```bash
# Use configuration in config.yaml
Rscript scripts/step3_cluster_seurat.R --config config.yaml

# Manually specify input file
Rscript scripts/step3_cluster_seurat.R \
  --config config.yaml \
  --input results/scp_analysis_preprocessed.tsv \
  --output_dir results/step3
```

#### Run Step 4 individually (Scanpy clustering)

```bash
# Use configuration in config.yaml
python scripts/step4_cluster_scanpy.py --config config.yaml

# Manually specify input file
python scripts/step4_cluster_scanpy.py \
  --config config.yaml \
  --input results/scp_analysis_preprocessed.tsv \
  --output_dir results/step4
```

Command-line parameters:

| Parameter | Description |
|-----------|-------------|
| `--config` | Configuration file path (default `config.yaml`) |
| `--input` | Input file path (overrides config.yaml) |
| `--output_dir` | Output directory (overrides default path) |
| `--exp_design` | Experimental design file (Step 2 only) |
| `--resolution` | Clustering resolution (Step 4 only, overrides config value) |


Sample column name format example:
```
[1] AraHypho_control_JAL_biaopi_1Cell_20251030_DIA_VA9000_5ng_Slot1-77_1_7896.d.PG.Quantity
```

## Output

```
results/
├── step1/
│   └── scp_analysis_preprocessed.tsv    # Preprocessed data
├── step2/
│   ├── 1se_plot.pdf                     # Protein identification statistics
│   ├── 2norm_plot.pdf                   # Normalization boxplot
│   ├── 3missval_plot.pdf                # Missing value visualization
│   ├── 4umap_plot.pdf                   # UMAP dimensionality reduction
│   ├── 5pca_plot.pdf                    # PCA principal components
│   ├── 6heatmap_volcano.pdf              # Heatmap + volcano plot
│   ├── 7data_results.csv                # Differential analysis results
│   ├── 6kmeanscluster.csv               # K-means clustering results
│   └── data.RData                       # R environment save
├── step3_seurat/
│   ├── 1protein_count.pdf               # Protein count per cell
│   ├── 2elbow_plot.pdf                  # PCA elbow plot
│   ├── 3umap.pdf                        # UMAP clustering
│   ├── 4markers.csv                     # Marker protein list
│   ├── 5top_markers.csv                 # Top marker proteins
│   ├── 6markers_dotplot_heatmap.pdf      # Dotplot + heatmap
│   ├── metadata.csv                     # Cell metadata
│   └── seurat_result.RData              # Seurat object
└── step4_scanpy/
    ├── umap.h5ad                        # Scanpy AnnData object
    ├── markers_all.csv                  # All marker proteins
    ├── metadata.csv                     # Cell metadata
    └── *.pdf                            # UMAP, dotplot, etc.
```

