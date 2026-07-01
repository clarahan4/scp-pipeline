#!/usr/bin/env python3
"""
Step 1: Spectronaut Single-cell proteomics data preprocessing

Functions:
1. Remove columns ending with SequencesIdentified
2. Remove specified metadata columns (PG.Qvalue, PG.Organisms, PG.ProteinDescriptions)
3. Remove protein rows starting with "Cont_"
4. Remove rows with protein or gene IDs containing special characters
5. Exclude sample columns matching specified patterns (e.g., 30Cell)
6. Remove protein rows where all sample values are NaN
7. Rename sample columns to concise format (e.g., control_1, cold_2)

Usage:
    python step1_preprocess.py --config config.yaml
    python step1_preprocess.py --config config.yaml --input raw_data.tsv --output processed.tsv
"""

import csv
import os
import re
import sys
import argparse
import yaml


def load_config(config_path):
    """Load YAML configuration file"""
    with open(config_path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


def has_special_characters(text):
    """Check if text contains comma, space, or quotes"""
    if text is None or text == '':
        return False
    return any(char in text for char in [',', ' ', '"', "'"])


def detect_sample_columns(columns):
    """Automatically detect sample columns (match Spectronaut PG.Quantity suffix)"""
    return [col for col in columns if col.endswith('.PG.Quantity')]


def drop_columns_by_config(columns, config):
    """Remove specified columns according to config, return kept columns"""
    drop_patterns = config.get('preprocess', {}).get('drop_columns', [])
    keep = []
    for col in columns:
        drop = False
        for pattern in drop_patterns:
            if pattern == "SequencesIdentified":
                if col.endswith("SequencesIdentified"):
                    drop = True
                    break
            elif col == pattern:
                drop = True
                break
        if not drop:
            keep.append(col)
    return keep


def rename_samples(sample_cols, config):
    """Rename sample columns according to regex pattern in config"""
    rename_pattern = config.get('preprocess', {}).get('rename_pattern', '')
    rename_template = config.get('preprocess', {}).get('rename_template', None)

    if not rename_pattern:
        return sample_cols

    new_names = []
    warn_count = 0
    for col in sample_cols:
        match = re.search(rename_pattern, col)
        if match and rename_template:
            try:
                # Pass all capture groups to format, numbered as group(1), group(2), ...
                # Also support named parameters: condition, number, cell_type
                groups = {}
                for i in range(len(match.groups()) + 1):
                    if i == 0:
                        continue
                    groups[f'group{i}'] = match.group(i)
                # Common aliases
                if len(match.groups()) >= 1:
                    groups['number'] = match.group(1)
                if len(match.groups()) >= 2:
                    groups['condition'] = match.group(2)
                if len(match.groups()) >= 3:
                    groups['cell_type'] = match.group(3)

                new_name = rename_template.format(**groups)
                new_names.append(new_name)
            except (IndexError, KeyError) as e:
                warn_count += 1
                if warn_count <= 3:
                    print(f"  Warning: Cannot rename '{col[:60]}...' ({e}), keeping original name")
                new_names.append(col)
        else:
            warn_count += 1
            if warn_count <= 3:
                print(f"  Warning: Column name '{col[:60]}...' does not match rename pattern, keeping original name")
            new_names.append(col)

    if warn_count > 3:
        print(f"  ... Total {warn_count} column names not matched")
    elif warn_count == 0:
        print(f"  All {len(sample_cols)} sample columns renamed successfully")

    return new_names


def main():
    parser = argparse.ArgumentParser(description='Step 1: Spectronaut Single-cell proteomics data preprocessing')
    parser.add_argument('--config', required=True, help='Configuration file path (YAML)')
    parser.add_argument('--input', default=None, help='Input file path (override config)')
    parser.add_argument('--output', default=None, help='Output file path (override config)')
    args = parser.parse_args()

    config = load_config(args.config)

    input_file = args.input or config.get('input', {}).get('raw_data')
    output_file = args.output or os.path.join(
        config.get('output', {}).get('dir', 'results'),
        config.get('output', {}).get('prefix', 'scp_analysis') + '_preprocessed.tsv'
    )

    if not input_file:
        print("Error: Input file not specified. Please specify via --input parameter or input.raw_data in config.yaml.")
        sys.exit(1)

    # Ensure output directory exists
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    print(f"{'='*60}")
    print(f"Step 1: Spectronaut Data preprocessing")
    print(f"{'='*60}")
    print(f"Input file: {input_file}")
    print(f"Output file: {output_file}")

    # ==========================================
    # Read data
    # ==========================================
    print(f"\n[1] Reading raw data...")
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found")
        sys.exit(1)

    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        rows = list(reader)

    print(f"  Original dimensions: {len(rows)} rows x {len(headers)} columns")

    # ==========================================
    # Detect sample columns
    # ==========================================
    sample_cols = detect_sample_columns(headers)
    meta_cols = [col for col in headers if col not in sample_cols]
    print(f"  Detected {len(sample_cols)} sample columns, {len(meta_cols)} metadata columns")

    # ==========================================
    # Remove specified columns
    # ==========================================
    print(f"\n[2] Removing specified columns...")
    keep_cols = drop_columns_by_config(headers, config)
    removed_count = len(headers) - len(keep_cols)
    if removed_count > 0:
        print(f"  Removed {removed_count} columns")
    else:
        print(f"  No columns to remove")

    # Rebuild data and column names
    col_index = {col: i for i, col in enumerate(headers)}
    headers = [h for h in keep_cols]
    keep_indices = [col_index[h] for h in keep_cols]
    rows = [[row[i] for i in keep_indices] for row in rows]

    # Re-detect sample columns
    sample_cols = detect_sample_columns(headers)
    meta_cols = [col for col in headers if col not in sample_cols]
    print(f"  Current dimensions: {len(rows)} rows x {len(headers)} columns (sample columns: {len(sample_cols)})")

    # Find key column indices
    protein_groups_idx = None
    genes_idx = None
    for i, h in enumerate(headers):
        if h == 'PG.ProteinGroups':
            protein_groups_idx = i
        elif h == 'PG.Genes':
            genes_idx = i

    # ==========================================
    # Remove proteins starting with Cont_
    # ==========================================
    if protein_groups_idx is not None:
        print(f"\n[3] Removing Cont_ proteins...")
        initial = len(rows)
        rows = [row for row in rows if not row[protein_groups_idx].startswith('Cont_')]
        print(f"  Removed {initial - len(rows)} Cont_ proteins")
    else:
        print(f"\n[3] PG.ProteinGroups column not found, skipping Cont_ removal")

    # ==========================================
    # Handle protein/gene ID rows with special characters
    # Strictly follow process_protein_data.py logic:
    # - Special characters = comma, space, quotes (semicolon not special, protein groups kept)
    # - Protein ID with special chars -> remove entire row
    # - Gene name with special chars -> clean and replace (space→_, quotes deleted) then keep
    # ==========================================
    print(f"\n[4] Processing rows with special characters...")
    initial = len(rows)
    filtered = []
    cleaned_count = 0
    for row in rows:
        skip = False
        # Protein ID contains comma/space/quotes -> remove entire row (semicolon not counted)
        if protein_groups_idx is not None and has_special_characters(row[protein_groups_idx]):
            skip = True
        if not skip:
            # Gene name contains special chars -> clean and replace then keep
            if genes_idx is not None and has_special_characters(row[genes_idx]):
                gene = row[genes_idx]
                gene = gene.replace(' ', '_').replace('"', '').replace("'", '')
                row[genes_idx] = gene
                cleaned_count += 1
            filtered.append(row)
    rows = filtered
    print(f"  Removed {initial - len(rows)} rows with special characters in protein ID")
    print(f"  Cleaned {cleaned_count} gene names (spaces/quotes removed, kept)")

    # ==========================================
    # Exclude sample columns matching specified patterns
    # ==========================================
    exclude_patterns = config.get('preprocess', {}).get('exclude_sample_patterns', [])
    if exclude_patterns:
        print(f"\n[5] Excluding sample columns matching patterns: {exclude_patterns}")
        remaining_sample_cols = []
        excluded_cols = []
        for col in sample_cols:
            excluded = False
            for pattern in exclude_patterns:
                if pattern in col:
                    excluded = True
                    break
            if excluded:
                excluded_cols.append(col)
            else:
                remaining_sample_cols.append(col)
        print(f"  Excluded {len(excluded_cols)} sample columns, kept {len(remaining_sample_cols)}")
    else:
        remaining_sample_cols = sample_cols
        print(f"\n[5] No sample columns to exclude (keeping all {len(remaining_sample_cols)})")

    # Also remove excluded columns from headers and delete corresponding column values from data rows
    excluded_set = set(excluded_cols) if exclude_patterns else set()
    if excluded_set:
        # Current indices for headers (with complete columns)
        col_index = {col: i for i, col in enumerate(headers)}
        keep_indices = [i for i, col in enumerate(headers) if col not in excluded_set]
        headers = [h for h in headers if h not in excluded_set]
        # Delete corresponding columns from each row
        rows = [[row[i] for i in keep_indices] for row in rows]

    # ==========================================
    # Remove rows where all sample values are NaN
    # ==========================================
    drop_all_nan = config.get('preprocess', {}).get('drop_all_nan_rows', True)
    if drop_all_nan:
        print(f"\n[6] Removing rows where all sample values are NaN...")
        remaining_sample_indices = [headers.index(col) for col in remaining_sample_cols]
        initial = len(rows)
        filtered = []
        for row in rows:
            all_nan = all(
                row[i] == 'NaN' or row[i] == 'nan' or row[i] == 'NA' or row[i] == ''
                for i in remaining_sample_indices
            )
            if not all_nan:
                filtered.append(row)
        rows = filtered
        print(f"  Removed {initial - len(rows)} all-NaN rows")
    else:
        print(f"\n[6] Skipping all-NaN row removal")

    print(f"  Current dimensions: {len(rows)} rows x {len(headers)} columns")

    # ==========================================
    # Rename sample columns
    # ==========================================
    print(f"\n[7] Renaming sample columns...")
    new_names = rename_samples(remaining_sample_cols, config)
    rename_map = dict(zip(remaining_sample_cols, new_names))
    headers = [rename_map.get(h, h) for h in headers]
    print(f"  Renamed {len(rename_map)} sample columns")

    # Examples
    if new_names:
        print(f"  Rename examples: {new_names[0]}, {new_names[1]}, {new_names[2]}")

    # ==========================================
    # Reorder columns: PG.ProteinGroups, PG.Genes, sample columns
    # ==========================================
    final_order = []
    for col in ['PG.ProteinGroups', 'PG.Genes']:
        if col in headers:
            final_order.append(col)
    for col in headers:
        if col not in final_order:
            final_order.append(col)

    col_index = {col: i for i, col in enumerate(headers)}
    final_indices = [col_index[c] for c in final_order]
    headers = final_order
    rows = [[row[i] for i in final_indices] for row in rows]

    # ==========================================
    # Save results
    # ==========================================
    print(f"\n[8] Saving results...")
    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
    with open(output_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(headers)
        writer.writerows(rows)

    print(f"  Output file: {output_file}")
    print(f"  Final dimensions: {len(rows)} rows x {len(headers)} columns")
    print(f"\nFinal column names:")
    for i, col in enumerate(headers):
        if i < 5 or i >= len(headers) - 2:
            print(f"  {i+1}. {col}")
        elif i == 5:
            print(f"  ... (omitting {len(headers) - 7} columns) ...")

    print(f"\n{'='*60}")
    print(f"Step 1 complete!")
    print(f"{'='*60}")

    return output_file


if __name__ == "__main__":
    main()
