#!/usr/bin/env python3
"""
Step 4: Single-cell proteomicsClustering analysis (Scanpy)

Use scanpy for single-cell proteomics dimensionality reduction, clustering, and differential analysis.

Usage:
    python step4_cluster_scanpy.py --config config.yaml
    python step4_cluster_scanpy.py --config config.yaml --input preprocessed.tsv --output_dir results/step4_scanpy
"""

import os
import sys
import argparse
import yaml
import pandas as pd
import numpy as np


def load_config(config_path):
    with open(config_path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser(description='Step 4: Scanpy Single-cell Clustering analysis')
    parser.add_argument('--config', required=True, help='Configuration file path (YAML)')
    parser.add_argument('--input', default=None, help='Preprocessed TSV file')
    parser.add_argument('--output_dir', default=None, help='Output directory')
    parser.add_argument('--resolution', type=float, default=None, help='Leiden clustering resolution (override config file)')
    args = parser.parse_args()

    config = load_config(args.config)

    # Prioritize command-line arguments, then use config, finally use default path
    if args.input:
        input_file = args.input
    elif config.get('scanpy', {}).get('input_file'):
        input_file = config['scanpy']['input_file']
    else:
        input_file = os.path.join(
            config.get('output', {}).get('dir', 'results'),
            config.get('output', {}).get('prefix', 'scp_analysis') + '_preprocessed.tsv'
        )
    output_dir = args.output_dir or os.path.join(
        config.get('output', {}).get('dir', 'results'), 'step4_scanpy'
    )
    scanpy_cfg = config.get('scanpy', {})

    if not os.path.exists(input_file):
        print(f"Error: Preprocessed file not found: {input_file}")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print("Step 4: Scanpy Single-cell Clustering analysis")
    print("=" * 60)
    print(f"Input data: {input_file}")
    print(f"Output directory: {output_dir}")
    print()

    import scanpy as sc
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt

    sc.settings.figdir = output_dir
    sc.settings.set_figure_params(dpi=150, figsize=(6, 5))

    # ==========================================
    # Reading data
    # ==========================================
    print("[1] Reading data...")
    protein_matrix = pd.read_csv(input_file, sep='\t', index_col=0)

    # Extract sample columns (exclude PG.Genes if present)
    meta_cols = [c for c in protein_matrix.columns if c in ('PG.Genes',)]
    sample_cols = [c for c in protein_matrix.columns if c not in meta_cols]

    # Create row names: ProteinID_GeneID
    if 'PG.Genes' in protein_matrix.columns:
        index_values = protein_matrix.index.astype(str) + '_' + protein_matrix['PG.Genes'].astype(str)
    else:
        index_values = protein_matrix.index.astype(str)

    expr = protein_matrix[sample_cols].copy()
    expr.index = index_values

    # Missing value imputation
    na_method = scanpy_cfg.get('na_fill_method', 'zero')
    if na_method == "zero":
        expr = expr.fillna(0)
        expr = expr.replace('NaN', 0).replace('nan', 0)
        print(f"  Missing value imputation: 0")
    elif na_method == "min":
        expr = expr.apply(lambda row: row.fillna(row.min()), axis=1)
        print(f"  Missing value imputation: row minimum values")

    print(f"  Protein matrix: {expr.shape[0]} proteins x {expr.shape[1]} cells")

    # Transpose: rows as cells, columns as proteins
    adata = sc.AnnData(expr.T)

    # ==========================================
    # Parse condition from column names
    # ==========================================
    print("\n[2] Parse sample information...")
    cell_names = adata.obs_names.tolist()
    conditions = []
    for name in cell_names:
        parts = name.split('_')
        conditions.append(parts[0] if parts else 'unknown')
    adata.obs['condition'] = pd.Categorical(conditions)

    print(f"  Condition statistics:")
    for cond, count in pd.Series(conditions).value_counts().items():
        print(f"    {cond}: {count}")

    # ==========================================
    # QC metrics calculation
    # ==========================================
    print("\n[3] QC metrics calculation...")
    sc.pp.calculate_qc_metrics(adata, qc_vars=[], percent_top=None,
                               log1p=False, inplace=True)

    print(f"  Proteins detected per cell: min={adata.obs['n_genes_by_counts'].min()}, "
          f"max={adata.obs['n_genes_by_counts'].max()}, "
          f"median={adata.obs['n_genes_by_counts'].median():.0f}")
    print(f"  Total proteins per cell: min={adata.obs['total_counts'].min():.1f}, "
          f"max={adata.obs['total_counts'].max():.1f}")

    # QC boxplot
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'],
                  jitter=0.4, multi_panel=True, save='_qc_violin.pdf')

    # ==========================================
    # Filtering
    # ==========================================
    print("\n[4] Filtering...")
    min_genes = scanpy_cfg.get('min_genes', 390)
    min_cells = scanpy_cfg.get('min_cells', 3)

    n_cells_before = adata.n_obs
    n_vars_before = adata.n_vars

    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    print(f"  cells: {n_cells_before} → {adata.n_obs} (Filtering {n_cells_before - adata.n_obs})")
    print(f"  proteins: {n_vars_before} → {adata.n_vars} (Filtering {n_vars_before - adata.n_vars})")

    # ==========================================
    # Normalization
    # ==========================================
    print("\n[5] Normalization...")
    target_sum = scanpy_cfg.get('target_sum', 10000)
    sc.pp.normalize_total(adata, target_sum=target_sum)
    print(f"  Normalize (target_sum={target_sum})")

    sc.pp.log1p(adata)
    print("  Log1p transform")

    # ==========================================
    # Highly variable gene selection
    # ==========================================
    print("\n[6] Highly variable gene selection...")
    hv_min_mean = scanpy_cfg.get('hv_min_mean', 0.0125)
    hv_max_mean = scanpy_cfg.get('hv_max_mean', 4.0)
    hv_min_disp = scanpy_cfg.get('hv_min_disp', 0.1)

    sc.pp.highly_variable_genes(adata, min_mean=hv_min_mean,
                                max_mean=hv_max_mean, min_disp=hv_min_disp)
    print(f"  Number of highly variable genes: {adata.var['highly_variable'].sum()}")

    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]

    # ==========================================
    # Scale, PCA, Neighbors
    # ==========================================
    print("\n[7] Scale, PCA, Neighbors...")
    scale_max = scanpy_cfg.get('scale_max_value', 10)
    sc.pp.scale(adata, max_value=scale_max)

    n_pcs = scanpy_cfg.get('n_pcs', 50)
    n_neighbors = scanpy_cfg.get('n_neighbors', 30)
    sc.tl.pca(adata, svd_solver="arpack", n_comps=min(n_pcs, adata.n_vars - 1))
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=min(n_pcs, adata.n_vars - 1))

    # PCA plot
    sc.pl.pca_variance_ratio(adata, log=True, save='_pca_variance.pdf')

    # ==========================================
    # UMAP
    # ==========================================
    print("\n[8] UMAP...")
    umap_min_dist = scanpy_cfg.get('umap_min_dist', 0.2)
    sc.tl.umap(adata, min_dist=umap_min_dist)

    sc.pl.umap(adata, color=['condition'], size=100,
               save='_by_condition.pdf')
    print("  Saved: umap_by_condition.pdf")

    # ==========================================
    # Leiden clustering
    # ==========================================
    print("\n[9] Leiden clustering...")
    resolution = args.resolution or scanpy_cfg.get('leiden_resolution', 0.5)
    print(f"  Use resolution: {resolution}")
    sc.tl.leiden(adata, resolution=resolution)

    n_clusters = len(adata.obs['leiden'].unique())
    print(f"  Identified {n_clusters} clusters (resolution={resolution})")

    sc.pl.umap(adata, color=['leiden'], size=100, legend_loc='on data',
               save='_clusters.pdf')
    print("  Saved: umap_clusters.pdf")

    # ==========================================
    # differentialproteinsanalysis
    # ==========================================
    print("\n[10] differentialproteinsanalysis...")
    marker_method = scanpy_cfg.get('marker_method', 'wilcoxon')
    top_n = scanpy_cfg.get('top_n_markers', 10)

    sc.tl.rank_genes_groups(adata, 'leiden', method=marker_method, pts=True)

    # Saved all marker results to CSV
    all_markers = sc.get.rank_genes_groups_df(adata, group=None)
    all_markers.to_csv(os.path.join(output_dir, 'markers_all.csv'), index=False)

    # Print top markers for each cluster
    for cluster in sorted(adata.obs["leiden"].unique()):
        markers = sc.get.rank_genes_groups_df(adata, group=cluster)
        top = markers.head(top_n)
        print(f"\n  Cluster {cluster} top {top_n} markers:")
        for _, row in top.iterrows():
            print(f"    {row['names']:20s}  logFC={row['logfoldchanges']:+.3f}  p={row['pvals_adj']:.2e}")

    # Dotplot and Violin
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=top_n, save='_dotplot.pdf')
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=top_n,
                                     values_to_plot='logfoldchanges',
                                     min_logfoldchange=2, vmax=7, vmin=-7,
                                     cmap='bwr', save='_dotplot_logfc.pdf')
    sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=top_n,
                                            cmap='viridis_r',
                                            save='_stacked_violin.pdf')
    print("  Saved: dotplot, stacked_violin PDF")

    # ==========================================
    # Savedresults
    # ==========================================
    print("\n[11] Savedresults...")
    adata.write_h5ad(os.path.join(output_dir, 'umap.h5ad'))

    # Saved metadata
    adata.obs.to_csv(os.path.join(output_dir, 'metadata.csv'))

    print(f"\n{'=' * 60}")
    print("Step 4 complete!")
    print(f"resultsSavedat: {output_dir}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
