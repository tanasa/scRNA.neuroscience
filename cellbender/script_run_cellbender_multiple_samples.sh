#!/bin/bash
# CellBender Sequential Pipeline for Temporal Lobe Samples
# GPU ONLY VERSION - Requires GPU/CUDA, exits if not available
# FIXED VERSION - Handles checkpoint errors gracefully
# Outputs saved to: /mnt/nfs/CX000008_DS1/projects/jaeyeon/h5_file_Dyselxia_r1/TEMPORAL_LOBE

# Configuration
CELLRANGER_OUTPUT_DIR="/mnt/nfs/CX000008_DS1/projects/jaeyeon/h5_file_Dyselxia_r1/TEMPORAL_LOBE"

# CellBender parameters
EPOCHS=100
MODEL="ambient"
LOW_COUNT_THRESHOLD=150
POSTERIOR_BATCH_SIZE=20
LEARNING_RATE=0.00005
USE_CUDA=true
CHECKPOINT_MINS=9999  # Delay checkpointing to avoid PyTorch pickle bug
RAM_THRESHOLD=90  # Stop if RAM usage exceeds this percentage

# Create output base directory if it doesn't exist
mkdir -p "$CELLRANGER_OUTPUT_DIR"

# Log file
LOG_FILE="${CELLRANGER_OUTPUT_DIR}/cellbender_pipeline_$(date +%Y%m%d_%H%M%S).log"

# Function to check RAM usage
check_ram() {
    local ram_usage=$(free | grep Mem | awk '{printf "%.0f", ($3/$2) * 100}')
    echo "$ram_usage"
}

# Function to display RAM status
show_ram_status() {
    local ram_usage=$(check_ram)
    local total_ram=$(free -g | grep Mem | awk '{print $2}')
    local used_ram=$(free -g | grep Mem | awk '{print $3}')
    local free_ram=$(free -g | grep Mem | awk '{print $4}')
    
    echo "RAM Status: ${ram_usage}% used (${used_ram}GB / ${total_ram}GB total, ${free_ram}GB free)" | tee -a "$LOG_FILE"
}

# Function to check CUDA availability
check_cuda() {
    if command -v nvidia-smi &> /dev/null; then
        echo "✓ nvidia-smi found" | tee -a "$LOG_FILE"
        nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader | tee -a "$LOG_FILE"
        
        # Check PyTorch CUDA
        python3 -c "import torch; print(f'PyTorch version: {torch.__version__}'); print(f'CUDA available: {torch.cuda.is_available()}'); print(f'CUDA version: {torch.version.cuda if torch.cuda.is_available() else \"N/A\"}'); print(f'GPU count: {torch.cuda.device_count() if torch.cuda.is_available() else 0}')" 2>&1 | tee -a "$LOG_FILE"
        
        if python3 -c "import torch; exit(0 if torch.cuda.is_available() else 1)" 2>/dev/null; then
            return 0
        else
            echo "⚠️  WARNING: CUDA not available in PyTorch, will use CPU" | tee -a "$LOG_FILE"
            return 1
        fi
    else
        echo "⚠️  WARNING: nvidia-smi not found, will use CPU" | tee -a "$LOG_FILE"
        return 1
    fi
}

echo "=========================================="
echo "CellBender Sequential Pipeline"
echo "Temporal Lobe Samples (16 samples)"
echo "GPU ONLY VERSION - Requires GPU/CUDA"
echo "FIXED VERSION - Handles checkpoint errors"
echo "=========================================="
echo "CellRanger output directory: $CELLRANGER_OUTPUT_DIR"
echo "CellBender parameters:"
echo "  Epochs: $EPOCHS"
echo "  Model: $MODEL"
echo "  Low count threshold: $LOW_COUNT_THRESHOLD"
echo "  Posterior batch size: $POSTERIOR_BATCH_SIZE"
echo "  Learning rate: $LEARNING_RATE"
echo "  CUDA requested: $USE_CUDA"
echo "  Checkpoint frequency: minimized (--checkpoint-mins $CHECKPOINT_MINS)"
echo "RAM threshold: ${RAM_THRESHOLD}% (will stop if exceeded)"
echo "Log file: $LOG_FILE"
echo "Start time: $(date)"
echo "=========================================="
echo "" | tee "$LOG_FILE"

# Check CUDA if requested - FORCE GPU ONLY
if [ "$USE_CUDA" = true ]; then
    echo "Checking CUDA availability..." | tee -a "$LOG_FILE"
    if check_cuda; then
        CUDA_FLAG="--cuda"
        echo "✓ Using CUDA (GPU acceleration)" | tee -a "$LOG_FILE"
    else
        echo "✗ ERROR: GPU/CUDA is required but not available!" | tee -a "$LOG_FILE"
        echo "✗ Please ensure:" | tee -a "$LOG_FILE"
        echo "  1. NVIDIA GPU is installed and drivers are loaded" | tee -a "$LOG_FILE"
        echo "  2. nvidia-smi command works" | tee -a "$LOG_FILE"
        echo "  3. PyTorch with CUDA support is installed" | tee -a "$LOG_FILE"
        echo "" | tee -a "$LOG_FILE"
        echo "Exiting..." | tee -a "$LOG_FILE"
        exit 1
    fi
    echo "" | tee -a "$LOG_FILE"
else
    CUDA_FLAG=""
    echo "⚠️  WARNING: USE_CUDA is set to false - using CPU (slower)" | tee -a "$LOG_FILE"
fi

# Initial RAM check
echo "Initial RAM check:" | tee -a "$LOG_FILE"
show_ram_status
echo "" | tee -a "$LOG_FILE"

# Define sample list
declare -a sample_list=(
    "ITG_4267_cellranger"
    "ITG_4285_cellranger"
    "ITG_4458_cellranger"
    "ITG_6145_cellranger"
    "MTG_4267_cellranger"
    "MTG_4285_cellranger"
    "MTG_4458_cellranger"
    "MTG_6145_cellranger"
    "STG_4267_cellranger"
    "STG_4285_cellranger"
    "STG_4458_cellranger"
    "STG_6145_cellranger"
    "tSVZ_4267_cellranger"
    "tSVZ_4285_cellranger"
    "tSVZ_6145_1_cellranger"
    "tSVZ_6145_2_cellranger"
)

echo "Sample list (16 samples):" | tee -a "$LOG_FILE"
for sample in "${sample_list[@]}"; do
    echo "  $sample" | tee -a "$LOG_FILE"
done
echo "" | tee -a "$LOG_FILE"

echo "Sample age groups:" | tee -a "$LOG_FILE"
echo "  35GW (35 gestational weeks): 4267 samples" | tee -a "$LOG_FILE"
echo "  7m (7 months postnatal): 4285, 4458 samples" | tee -a "$LOG_FILE"
echo "  39GW (39 gestational weeks): 6145 samples" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Function to run CellBender for a single sample
run_cellbender() {
    local sample_name=$1
    local cellranger_dir="${CELLRANGER_OUTPUT_DIR}/${sample_name}"
    local input_h5="${cellranger_dir}/outs/raw_feature_bc_matrix.h5"
    local output_h5="${cellranger_dir}/outs/raw_feature_bc_matrix.cellbender.h5"
    
    echo "" | tee -a "$LOG_FILE"
    echo "==========================================" | tee -a "$LOG_FILE"
    echo "Processing: $sample_name" | tee -a "$LOG_FILE"
    echo "CellRanger directory: $cellranger_dir" | tee -a "$LOG_FILE"
    echo "Input h5: $input_h5" | tee -a "$LOG_FILE"
    echo "Output h5: $output_h5" | tee -a "$LOG_FILE"
    echo "==========================================" | tee -a "$LOG_FILE"
    
    # Check if CellRanger directory exists
    if [ ! -d "$cellranger_dir" ]; then
        echo "✗ ERROR: CellRanger directory does not exist: $cellranger_dir" | tee -a "$LOG_FILE"
        return 1
    fi
    
    # Check if input file exists
    if [ ! -f "$input_h5" ]; then
        echo "✗ ERROR: Input h5 file does not exist: $input_h5" | tee -a "$LOG_FILE"
        return 1
    fi
    
    # Check if already completed
    if [ -f "$output_h5" ]; then
        echo "✓ Output file already exists. Skipping." | tee -a "$LOG_FILE"
        return 0
    fi
    
    # Show input file info
    echo "Input file info:" | tee -a "$LOG_FILE"
    ls -lh "$input_h5" | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"
    
    # Remove any existing checkpoint files to avoid hash mismatch
    local checkpoint_dir="${cellranger_dir}/cellbender_checkpoint"
    local checkpoint_tarball="${cellranger_dir}/ckpt.tar.gz"
    
    if [ -d "$checkpoint_dir" ]; then
        echo "Removing existing checkpoint directory: $checkpoint_dir" | tee -a "$LOG_FILE"
        rm -rf "$checkpoint_dir"
    fi
    
    if [ -f "$checkpoint_tarball" ]; then
        echo "Removing existing checkpoint tarball: $checkpoint_tarball" | tee -a "$LOG_FILE"
        rm -f "$checkpoint_tarball"
    fi
    
    # Log start
    start_time=$(date +%s)
    echo "Start time: $(date)" | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"
    
    # Build CellBender command
    echo "Running CellBender remove-background..." | tee -a "$LOG_FILE"
    if [ -n "$CUDA_FLAG" ]; then
        echo "Using CUDA (GPU acceleration)" | tee -a "$LOG_FILE"
    else
        echo "Using CPU" | tee -a "$LOG_FILE"
    fi
    
    # Construct command
    local cmd="cellbender remove-background"
    cmd="$cmd --input \"$input_h5\""
    cmd="$cmd --output \"$output_h5\""
    cmd="$cmd --epochs $EPOCHS"
    cmd="$cmd --model $MODEL"
    cmd="$cmd --low-count-threshold $LOW_COUNT_THRESHOLD"
    cmd="$cmd --posterior-batch-size $POSTERIOR_BATCH_SIZE"
    cmd="$cmd --learning-rate $LEARNING_RATE"
    cmd="$cmd --checkpoint-mins $CHECKPOINT_MINS"
    if [ -n "$CUDA_FLAG" ]; then
        cmd="$cmd $CUDA_FLAG"
    fi
    
    echo "Command: $cmd" | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"
    
    # Run CellBender - capture output but don't exit on error
    # The checkpoint saving may fail, but output file might still be created
    eval "$cmd" 2>&1 | tee -a "$LOG_FILE"
    local cellbender_exit_code=${PIPESTATUS[0]}
    
    # Check if output file was created (even if CellBender exited with error)
    if [ -f "$output_h5" ]; then
        echo "" | tee -a "$LOG_FILE"
        echo "✓ SUCCESS: Output file created: $output_h5" | tee -a "$LOG_FILE"
        
        # Check file size to ensure it's valid
        local file_size=$(stat -f%z "$output_h5" 2>/dev/null || stat -c%s "$output_h5" 2>/dev/null)
        if [ "$file_size" -gt 1000 ]; then
            echo "✓ Output file size: $(numfmt --to=iec-i --suffix=B $file_size 2>/dev/null || echo "${file_size} bytes")" | tee -a "$LOG_FILE"
            end_time=$(date +%s)
            duration=$((end_time - start_time))
            hours=$((duration / 3600))
            minutes=$(((duration % 3600) / 60))
            echo "✓ Completed in ${hours}h ${minutes}m" | tee -a "$LOG_FILE"
            echo "End time: $(date)" | tee -a "$LOG_FILE"
            return 0
        else
            echo "⚠️  WARNING: Output file exists but is suspiciously small" | tee -a "$LOG_FILE"
            return 1
        fi
    else
        echo "" | tee -a "$LOG_FILE"
        echo "✗ FAILED: Output file was not created" | tee -a "$LOG_FILE"
        echo "CellBender exit code: $cellbender_exit_code" | tee -a "$LOG_FILE"
        end_time=$(date +%s)
        duration=$((end_time - start_time))
        hours=$((duration / 3600))
        minutes=$(((duration % 3600) / 60))
        echo "Duration: ${hours}h ${minutes}m" | tee -a "$LOG_FILE"
        return 1
    fi
}

# Verify CellRanger output directory exists
if [ ! -d "$CELLRANGER_OUTPUT_DIR" ]; then
    echo "✗ ERROR: CellRanger output directory does not exist: $CELLRANGER_OUTPUT_DIR" | tee -a "$LOG_FILE"
    exit 1
fi

echo "Verifying CellRanger sample directories..." | tee -a "$LOG_FILE"
missing_samples=0
for sample in "${sample_list[@]}"; do
    if [ ! -d "${CELLRANGER_OUTPUT_DIR}/${sample}" ]; then
        echo "⚠️  WARNING: Directory not found: ${CELLRANGER_OUTPUT_DIR}/${sample}" | tee -a "$LOG_FILE"
        ((missing_samples++))
    fi
done

if [ $missing_samples -gt 0 ]; then
    echo "" | tee -a "$LOG_FILE"
    echo "⚠️  WARNING: $missing_samples sample directories not found" | tee -a "$LOG_FILE"
    echo "Continuing with available samples..." | tee -a "$LOG_FILE"
fi
echo "" | tee -a "$LOG_FILE"

# Run each sample sequentially
total_samples=${#sample_list[@]}
completed=0
failed=0
skipped=0

pipeline_start=$(date +%s)

for i in "${!sample_list[@]}"; do
    sample_num=$((i + 1))
    sample="${sample_list[$i]}"
    
    # === RAM CHECK BEFORE STARTING EACH SAMPLE ===
    echo "" | tee -a "$LOG_FILE"
    echo "=========================================="  | tee -a "$LOG_FILE"
    echo "PRE-SAMPLE RAM CHECK"  | tee -a "$LOG_FILE"
    echo "=========================================="  | tee -a "$LOG_FILE"
    
    current_ram=$(check_ram)
    show_ram_status
    
    if [ "$current_ram" -gt "$RAM_THRESHOLD" ]; then
        echo "" | tee -a "$LOG_FILE"
        echo "🛑 CRITICAL: RAM usage (${current_ram}%) exceeds threshold (${RAM_THRESHOLD}%)" | tee -a "$LOG_FILE"
        echo "🛑 STOPPING PIPELINE FOR SAFETY" | tee -a "$LOG_FILE"
        echo "" | tee -a "$LOG_FILE"
        echo "Completed samples so far: $completed" | tee -a "$LOG_FILE"
        echo "Remaining samples: $((total_samples - completed - failed - skipped))" | tee -a "$LOG_FILE"
        echo "" | tee -a "$LOG_FILE"
        echo "To resume:" | tee -a "$LOG_FILE"
        echo "1. Wait for RAM to be released (check with 'free -g')" | tee -a "$LOG_FILE"
        echo "2. Rerun this script - it will skip completed samples" | tee -a "$LOG_FILE"
        echo "" | tee -a "$LOG_FILE"
        
        # Calculate elapsed time
        pipeline_end=$(date +%s)
        total_duration=$((pipeline_end - pipeline_start))
        total_hours=$((total_duration / 3600))
        total_minutes=$(((total_duration % 3600) / 60))
        
        echo "Elapsed time: ${total_hours}h ${total_minutes}m" | tee -a "$LOG_FILE"
        echo "End time: $(date)" | tee -a "$LOG_FILE"
        echo "" | tee -a "$LOG_FILE"
        
        exit 2  # Exit code 2 = stopped due to high RAM
    fi
    
    echo "" | tee -a "$LOG_FILE"
    echo "=========================================="  | tee -a "$LOG_FILE"
    echo "SAMPLE $sample_num of $total_samples: $sample"  | tee -a "$LOG_FILE"
    echo "=========================================="  | tee -a "$LOG_FILE"
    
    if [ ! -d "${CELLRANGER_OUTPUT_DIR}/${sample}" ]; then
        echo "⚠️  CellRanger directory not found. Skipping..." | tee -a "$LOG_FILE"
        ((skipped++))
        continue
    fi
    
    if run_cellbender "$sample"; then
        ((completed++))
    else
        ((failed++))
        echo "" | tee -a "$LOG_FILE"
        echo "⚠️  Sample failed. Continuing with next sample..." | tee -a "$LOG_FILE"
    fi
    
    echo "" | tee -a "$LOG_FILE"
    echo "Progress: $completed completed, $failed failed, $skipped skipped, $((total_samples - completed - failed - skipped)) remaining" | tee -a "$LOG_FILE"
    
    # Post-sample RAM check
    echo "" | tee -a "$LOG_FILE"
    echo "Post-sample RAM check:" | tee -a "$LOG_FILE"
    show_ram_status
    
    # Wait between samples to ensure clean RAM release
    if [ $sample_num -lt $total_samples ]; then
        echo "" | tee -a "$LOG_FILE"
        echo "Waiting 30 seconds for RAM cleanup before next sample..." | tee -a "$LOG_FILE"
        sleep 30
    fi
done

# Calculate total time
pipeline_end=$(date +%s)
total_duration=$((pipeline_end - pipeline_start))
total_hours=$((total_duration / 3600))
total_minutes=$(((total_duration % 3600) / 60))

# Final RAM check
echo "" | tee -a "$LOG_FILE"
echo "Final RAM check:" | tee -a "$LOG_FILE"
show_ram_status

# Final summary
echo "" | tee -a "$LOG_FILE"
echo "==========================================" | tee -a "$LOG_FILE"
echo "PIPELINE COMPLETE" | tee -a "$LOG_FILE"
echo "==========================================" | tee -a "$LOG_FILE"
echo "Total samples: $total_samples" | tee -a "$LOG_FILE"
echo "Completed: $completed" | tee -a "$LOG_FILE"
echo "Failed: $failed" | tee -a "$LOG_FILE"
echo "Skipped: $skipped" | tee -a "$LOG_FILE"
echo "Total time: ${total_hours}h ${total_minutes}m" | tee -a "$LOG_FILE"
echo "End time: $(date)" | tee -a "$LOG_FILE"
echo "Output directory: $CELLRANGER_OUTPUT_DIR" | tee -a "$LOG_FILE"
echo "==========================================" | tee -a "$LOG_FILE"

# List completed samples
if [ $completed -gt 0 ]; then
    echo "" | tee -a "$LOG_FILE"
    echo "Completed samples:" | tee -a "$LOG_FILE"
    for sample in "${sample_list[@]}"; do
        output_h5="${CELLRANGER_OUTPUT_DIR}/${sample}/outs/raw_feature_bc_matrix.cellbender.h5"
        if [ -f "$output_h5" ]; then
            echo "  ✓ $sample" | tee -a "$LOG_FILE"
            echo "    → $output_h5" | tee -a "$LOG_FILE"
        fi
    done
fi

# List failed samples
if [ $failed -gt 0 ]; then
    echo "" | tee -a "$LOG_FILE"
    echo "Failed samples:" | tee -a "$LOG_FILE"
    for sample in "${sample_list[@]}"; do
        output_h5="${CELLRANGER_OUTPUT_DIR}/${sample}/outs/raw_feature_bc_matrix.cellbender.h5"
        if [ ! -f "$output_h5" ]; then
            echo "  ✗ $sample" | tee -a "$LOG_FILE"
        fi
    done
fi

echo "" | tee -a "$LOG_FILE"
echo "Log file saved to: $LOG_FILE" | tee -a "$LOG_FILE"

if [ $failed -gt 0 ]; then
    echo "" | tee -a "$LOG_FILE"
    echo "⚠️  Some samples failed. Check individual logs for details." | tee -a "$LOG_FILE"
    exit 1
else
    echo "" | tee -a "$LOG_FILE"
    echo "✓ All samples completed successfully!" | tee -a "$LOG_FILE"
    exit 0
fi

