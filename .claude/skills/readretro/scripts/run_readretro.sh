#!/bin/bash
# READRetro Helper Script
# Usage: bash run_readretro.sh "SMILES_STRING" [options]

set -e

# Default parameters
SMILES=""
GPU_ID=0
MODE="single"
OUTPUT_DIR="result"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -s|--smiles)
            SMILES="$2"
            shift 2
            ;;
        -g|--gpu)
            GPU_ID="$2"
            shift 2
            ;;
        -m|--mode)
            MODE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: bash run_readretro.sh [options]"
            echo ""
            echo "Options:"
            echo "  -s, --smiles SMILES    SMILES string of target molecule (required for single mode)"
            echo "  -g, --gpu GPU_ID       GPU device ID (default: 0)"
            echo "  -m, --mode MODE        Mode: single|batch|eval (default: single)"
            echo "  -o, --output DIR       Output directory (default: result)"
            echo "  -h, --help            Show this help message"
            echo ""
            echo "Examples:"
            echo "  bash run_readretro.sh -s 'O=C1C=C2C=CC(O)CC2O1'"
            echo "  bash run_readretro.sh -m batch -g 0"
            echo "  bash run_readretro.sh -m eval"
            exit 0
            ;;
        *)
            SMILES="$1"
            shift
            ;;
    esac
done

# Validate inputs
if [[ "$MODE" == "single" && -z "$SMILES" ]]; then
    echo "Error: SMILES string required for single mode"
    echo "Usage: bash run_readretro.sh -s 'SMILES_STRING'"
    exit 1
fi

# Activate conda environment
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate readretro

# Change to READRetro directory
cd /home/xiejinxiang/chemoenzymatic-retrosynthesis-agent/READRetro

# Run based on mode
case $MODE in
    single)
        echo "Running READRetro for molecule: $SMILES"
        CUDA_VISIBLE_DEVICES=$GPU_ID python run.py "$SMILES"
        ;;
    batch)
        echo "Running READRetro in batch mode..."
        CUDA_VISIBLE_DEVICES=$GPU_ID python run_mp.py
        ;;
    eval)
        echo "Running single-step evaluation..."
        CUDA_VISIBLE_DEVICES=$GPU_ID python eval_single.py
        ;;
    *)
        echo "Error: Unknown mode '$MODE'"
        echo "Valid modes: single, batch, eval"
        exit 1
        ;;
esac

echo "Done!"
