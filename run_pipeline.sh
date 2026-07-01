#!/usr/bin/env bash
# ============================================================
# Single-Cell Proteomics Pipeline - Sequential Runner
# ============================================================
# Usage:
#   bash run_pipeline.sh [config.yaml]
#   bash run_pipeline.sh [config.yaml] --step1 --step2 --step3
#   bash run_pipeline.sh [config.yaml] --step3 --step4   # only clustering
# ============================================================
set -euo pipefail

CONFIG="${1:-config.yaml}"

if [ ! -f "$CONFIG" ]; then
    echo "Error: Config file '$CONFIG' not found"
    echo "Usage: bash run_pipeline.sh [config.yaml]"
    exit 1
fi

echo "=============================================="
echo "  Single-Cell Proteomics Pipeline"
echo "=============================================="
echo "Config: $CONFIG"

# Read output dir
OUTPUT_DIR=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['output']['dir'])" 2>/dev/null || echo "results")

# Read step toggles from config (default: step1=true, step2/3/4 from config)
read_config() {
    python3 -c "
import yaml
c = yaml.safe_load(open('$CONFIG'))
steps = c.get('steps', {})
print(steps.get('step1', True))
print(steps.get('step2', True))
print(steps.get('step3', c.get('seurat', {}).get('run', True)))
print(steps.get('step4', c.get('scanpy', {}).get('run', False)))
"
}

# Command-line args override config
RUN_STEP1=true
RUN_STEP2=true
RUN_STEP3=true
RUN_STEP4=true

if [ $# -gt 1 ]; then
    RUN_STEP1=false; RUN_STEP2=false; RUN_STEP3=false; RUN_STEP4=false
    for arg in "${@:2}"; do
        case "$arg" in
            --step1) RUN_STEP1=True ;;
            --step2) RUN_STEP2=True ;;
            --step3) RUN_STEP3=True ;;
            --step4) RUN_STEP4=True ;;
        esac
    done
else
    read_config_vals=$(read_config)
    RUN_STEP1=$(echo "$read_config_vals" | sed -n '1p')
    RUN_STEP2=$(echo "$read_config_vals" | sed -n '2p')
    RUN_STEP3=$(echo "$read_config_vals" | sed -n '3p')
    RUN_STEP4=$(echo "$read_config_vals" | sed -n '4p')
fi

echo ""
echo "Steps to run:"
[ "$RUN_STEP1" = "True" ] && echo "  [ON]  Step 1: Data preprocessing" || echo "  [OFF] Step 1: Data preprocessing"
[ "$RUN_STEP2" = "True" ] && echo "  [ON]  Step 2: Differential expression (DEP)" || echo "  [OFF] Step 2: Differential expression (DEP)"
[ "$RUN_STEP3" = "True" ] && echo "  [ON]  Step 3: Seurat clustering" || echo "  [OFF] Step 3: Seurat clustering"
[ "$RUN_STEP4" = "True" ] && echo "  [ON]  Step 4: Scanpy clustering" || echo "  [OFF] Step 4: Scanpy clustering"
echo ""
echo "Output: $OUTPUT_DIR/"
echo ""

# Check input file
INPUT_FILE=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['input']['raw_data'])" 2>/dev/null)
if [ "$RUN_STEP1" = "True" ] && [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found"
    exit 1
fi

# ============================================
# Step 1: Preprocessing
# ============================================
if [ "$RUN_STEP1" = "True" ]; then
    echo ">>> [Step 1] Data preprocessing..."
    python3 scripts/step1_preprocess.py --config "$CONFIG"

    echo ">>> Generating experimental design..."
    python3 scripts/generate_exp_design.py --config "$CONFIG"
else
    # Check if preprocessed file exists (needed by downstream steps)
    PREPROC=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['output']['prefix'])" 2>/dev/null || echo "scp_analysis")
    if [ ! -f "$OUTPUT_DIR/${PREPROC}_preprocessed.tsv" ]; then
        echo "Error: Step 1 is disabled but preprocessed file not found:"
        echo "  $OUTPUT_DIR/${PREPROC}_preprocessed.tsv"
        echo "Either enable step1 or run it first."
        exit 1
    fi
    echo ">>> [Step 1] Skipped (using existing preprocessed data)"
fi

# ============================================
# Step 2: Differential Expression
# ============================================
if [ "$RUN_STEP2" = "True" ]; then
    echo ""
    echo ">>> [Step 2] Differential expression analysis..."
    mkdir -p "$OUTPUT_DIR/step2"
    Rscript scripts/step2_de_analysis.R --config "$CONFIG"
else
    echo ">>> [Step 2] Skipped"
fi

# ============================================
# Step 3: Seurat Clustering
# ============================================
if [ "$RUN_STEP3" = "True" ]; then
    echo ""
    echo ">>> [Step 3] Seurat clustering..."
    mkdir -p "$OUTPUT_DIR/step3_seurat"
    Rscript scripts/step3_cluster_seurat.R --config "$CONFIG"
else
    echo ">>> [Step 3] Skipped"
fi

# ============================================
# Step 4: Scanpy Clustering
# ============================================
if [ "$RUN_STEP4" = "True" ]; then
    echo ""
    echo ">>> [Step 4] Scanpy clustering..."
    mkdir -p "$OUTPUT_DIR/step4_scanpy"
    python3 scripts/step4_cluster_scanpy.py --config "$CONFIG"
else
    echo ">>> [Step 4] Skipped"
fi

echo ""
echo "=============================================="
echo "  Pipeline complete!"
echo "=============================================="
echo ""
