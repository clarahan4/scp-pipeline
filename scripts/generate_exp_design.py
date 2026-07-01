#!/usr/bin/env python3
"""
Automatically generate ExpDesign.csv experimental design file

Read sample column names from preprocessed TSV file,
extract condition information using regex pattern from config,
generate experimental design file required by DEP package.

Usage:
    python generate_exp_design.py --config config.yaml
    python generate_exp_design.py --input preprocessed.tsv --output exp_design.csv
"""

import csv
import os
import re
import argparse
import yaml


def load_config(config_path):
    with open(config_path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser(description='Automatically generate experimental design file ExpDesign.csv')
    parser.add_argument('--config', required=True, help='Configuration file path')
    parser.add_argument('--input', default=None, help='Preprocessed TSV file (override config)')
    parser.add_argument('--output', default=None, help='Output CSV file path (override config)')
    args = parser.parse_args()

    config = load_config(args.config)

    # Determine input file: if config exists, try to get preprocessed path from config
    input_file = args.input
    if not input_file:
        output_dir = config.get('output', {}).get('dir', 'results')
        prefix = config.get('output', {}).get('prefix', 'scp_analysis')
        input_file = os.path.join(output_dir, prefix + '_preprocessed.tsv')

    output_file = args.output or config.get('input', {}).get('exp_design', 'test/test_exp_design.csv')

    if not os.path.exists(input_file):
        print(f"Error: Preprocessed file '{input_file}' not found, please run step1 first.")
        return

    print(f"{'='*60}")
    print(f"Generate experimental design file")
    print(f"{'='*60}")
    print(f"Input file: {input_file}")
    print(f"Output file: {output_file}")

    # Read column names
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)

    # Exclude metadata columns (PG.ProteinGroups, PG.Genes), only keep sample columns
    sample_cols = [h for h in headers if h not in ('PG.ProteinGroups', 'PG.Genes')]
    print(f"Detected {len(sample_cols)} sample columns")

    # Get rename pattern from config to parse condition
    rename_pattern = config.get('preprocess', {}).get('rename_pattern', '')

    # If rename pattern exists, use same pattern to extract condition
    if rename_pattern:
        condition_map = {}
        for col in sample_cols:
            match = re.search(rename_pattern, col)
            # Try to extract from new format column names (e.g., control_1 → condition=control, number=1)
            # Also try to extract from original format column names
            if match:
                condition = match.group(1)  # condition
                condition_map[col] = condition
            else:
                # Try to parse from simplified name (e.g., control_1 → control)
                parts = col.rsplit('_', 1)
                if len(parts) == 2:
                    condition_map[col] = parts[0]
                else:
                    condition_map[col] = 'unknown'
    else:
        # Without pattern, use part before first underscore in column name as condition
        condition_map = {}
        for col in sample_cols:
            parts = col.split('_')
            condition_map[col] = parts[0] if parts else 'unknown'

    # Group by condition, assign replicate numbers
    condition_replicates = {}
    for col in sample_cols:
        cond = condition_map[col]
        if cond not in condition_replicates:
            condition_replicates[cond] = 0
        condition_replicates[cond] += 1

    # Generate rows
    rows = []
    condition_counters = {c: 0 for c in condition_replicates}
    for col in sample_cols:
        cond = condition_map[col]
        condition_counters[cond] += 1
        rows.append({
            'label': col,
            'condition': cond,
            'replicate': condition_counters[cond]
        })

    # Save
    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
    with open(output_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['label', 'condition', 'replicate'])
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nExperimental design summary:")
    for cond, count in sorted(condition_replicates.items()):
        print(f"  {cond}: {count} samples")

    print(f"\nOutput: {output_file} ({len(rows)} rows)")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
