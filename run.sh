#!/bin/bash
#SBATCH --job-name=multinichenet_analysis
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vishoza@uab.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=350G
#SBATCH --time=24:00:00
#SBATCH --partition=largemem
#SBATCH --output=multinichenet_analysis_%A_%a.out
#SBATCH --error=multinichenet_analysis_%A_%a.err

# Load necessary modules
module load Java/19.0.2
module load Anaconda3

# Set JAVA_HOME to the location of the loaded Java
export JAVA_HOME=$(dirname $(dirname $(which java)))

# Check JAVA_HOME environment variable
echo $JAVA_HOME

# Set up Nextflow
NXF_DIR="/data/user/vishoza/nextflow"
export PATH=$PATH:$NXF_DIR
export NXF_HOME=$HOME/.nextflow

# Set R memory size for larger datasets
export R_MAX_VSIZE=300000000000 # Set to 300GB

# Specify input and output directories
INPUT_DIR="/data/project/zindl_lab/pragun/garrett"
OUTPUT_DIR="/data/project/zindl_lab/H2Ab1/data/multinichenet_analysis/"

# Ensure the output directory exists
mkdir -p $OUTPUT_DIR

# Create subdirectories for Nextflow work and temporary files within OUTPUT_DIR
WORK_DIR="${OUTPUT_DIR}/work"
TMP_DIR="${OUTPUT_DIR}/tmp"
mkdir -p $OUTPUT_DIR $WORK_DIR $TMP_DIR

# Set environment variables for temporary directories
export TMPDIR=$TMP_DIR
export NXF_TEMP=$TMP_DIR

# Define input parameters
SEURAT_FILE="${INPUT_DIR}/seurat_object.rds"  # Path to your Seurat RDS file
CELL_TYPE_COLUMN="cell_type_subsets"  # Column in Seurat metadata containing cell types

# Force pull the latest version of the pipeline
echo "Pulling latest pipeline version..."
nextflow pull vishaloza/ccc-pipeline -r main

# Run the pipeline with MultiNicheNet
nextflow run vishaloza/ccc-pipeline \
  -r main \
  -latest \
  -resume \
  -profile conda,slurm \
  -work-dir $WORK_DIR \
  --input $SEURAT_FILE \
  --input_type rds_multi \
  --celltype_column $CELL_TYPE_COLUMN \
  --outdir $OUTPUT_DIR \
  -with-report \
  -with-trace \
  -with-timeline

echo "MultiNicheNet analysis complete!"