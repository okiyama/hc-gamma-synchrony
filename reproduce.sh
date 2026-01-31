#!/bin/bash
# reproduce.sh - Full reproduction of hippocampal-cortical gamma synchrony analysis
#
# This script will:
# 1. Create a fresh Python environment
# 2. Install all dependencies
# 3. Run the full analysis pipeline (~10-12 hours)
# 4. Generate all figures
#
# Usage: ./reproduce.sh
#
# Requirements:
# - Python 3.8+ or conda
# - ~50 GB disk space for Allen Institute data
# - ~16 GB RAM

set -e  # Exit on error

echo "============================================================"
echo "Hippocampal-Cortical Gamma Synchrony Analysis"
echo "Full Reproduction Script"
echo "============================================================"
echo ""

# Check Python version
PYTHON_VERSION=$(python3 --version 2>&1 | cut -d' ' -f2 | cut -d'.' -f1,2)
echo "Python version: $PYTHON_VERSION"

# Create virtual environment if it doesn't exist
if [ ! -d ".venv" ]; then
    echo ""
    echo "Creating virtual environment..."
    python3 -m venv .venv
fi

# Activate environment
echo "Activating environment..."
source .venv/bin/activate

# Install dependencies
echo ""
echo "Installing dependencies..."
pip install --upgrade pip
pip install -r requirements.txt

# Create output directories
mkdir -p results
mkdir -p figures

# Run analysis
echo ""
echo "============================================================"
echo "Starting analysis pipeline"
echo "This will take approximately 10-12 hours"
echo "============================================================"
echo ""

START_TIME=$(date +%s)

python analysis/run_analysis.py \
    --output-dir results \
    --n-surrogates 5000 \
    --seed 7829

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
HOURS=$((ELAPSED / 3600))
MINUTES=$(((ELAPSED % 3600) / 60))

echo ""
echo "Analysis complete in ${HOURS}h ${MINUTES}m"

# Generate figures
echo ""
echo "Generating figures..."
python figures/generate_figures.py --input results/all_results.json --output figures/

echo ""
echo "============================================================"
echo "REPRODUCTION COMPLETE"
echo "============================================================"
echo ""
echo "Results: results/all_results.json"
echo "Figures: figures/"
echo ""
echo "To verify results match published data:"
echo "  python analysis/verify_results.py"
echo ""
