!/bin/bash

# Simple CellBender Script
# Usage: ./run_cellbender_simple.sh <input_file>

set -e  # Exit on any error

# Check if input file is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <input_file>"
    echo "Example: $0 raw_feature_bc_matrix.h5"
    exit 1
fi

INPUT_FILE="$1"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' does not exist"
    exit 1
fi

# Get the directory and basename of the input file
INPUT_DIR=$(dirname "$INPUT_FILE")
INPUT_BASENAME=$(basename "$INPUT_FILE")
INPUT_NAME="${INPUT_BASENAME%.*}"  # Remove extension

# Create output filename
OUTPUT_FILE="${INPUT_DIR}/${INPUT_NAME}.cellbender.h5"

echo "Input:  $INPUT_FILE"
echo "Output: $OUTPUT_FILE"
echo "Running CellBender..."

# Run CellBender
cellbender remove-background \
    --input "$INPUT_FILE" \
    --output "$OUTPUT_FILE" \
    --epochs 100 \
    --model ambient \
    --low-count-threshold 150 \
    --posterior-batch-size 20 \
    --cuda \
    --learning-rate 0.00005

echo "Done! Output: $OUTPUT_FILE"




