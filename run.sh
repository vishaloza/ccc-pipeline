#!/bin/bash
# Helper script to run the LIANA + NicheNet pipeline with Singularity 3.5.2

# Default values
INPUT=""
SENDER=""
RECEIVER=""
OUTDIR="results"
PROFILE="singularity_3_5,slurm"  # Use the Singularity 3.5 compatible profile by default

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --input)
      INPUT="$2"
      shift 2
      ;;
    --sender)
      SENDER="$2"
      shift 2
      ;;
    --receiver)
      RECEIVER="$2"
      shift 2
      ;;
    --outdir)
      OUTDIR="$2"
      shift 2
      ;;
    --profile)
      PROFILE="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1"
      exit 1
      ;;
  esac
done

# Check for required arguments
if [[ -z "$INPUT" || -z "$SENDER" || -z "$RECEIVER" ]]; then
  echo "Error: Required arguments missing"
  echo "Usage: ./run_hpc.sh --input <h5ad_file> --sender <sender_celltype> --receiver <receiver_celltype> [--outdir <output_dir>] [--profile <nextflow_profile>]"
  exit 1
fi

# Check if file exists
if [[ ! -f "$INPUT" ]]; then
  echo "Error: Input file not found: $INPUT"
  exit 1
fi

# Check Singularity version
SINGULARITY_VERSION=$(singularity --version 2>/dev/null | grep -oP '(?<=version\s)[0-9]+\.[0-9]+(\.[0-9]+)?')
echo "Detected Singularity version: $SINGULARITY_VERSION"

# Set appropriate profile based on Singularity version
if [[ "$PROFILE" == *"singularity"* && "$SINGULARITY_VERSION" == "3.5"* ]]; then
  # Replace standard singularity profile with singularity_3_5
  PROFILE=$(echo "$PROFILE" | sed 's/singularity/singularity_3_5/g')
  echo "Using compatibility profile for Singularity 3.5.x: $PROFILE"
fi

# Create the work directory with appropriate permissions
mkdir -p "$OUTDIR"
mkdir -p work
chmod -R 777 work

# Export Singularity environment variables
export SINGULARITY_BINDPATH="$PWD:/data,$PWD:/workdir,$PWD:$PWD"
export SINGULARITY_CACHEDIR="$PWD/.singularity"
mkdir -p $SINGULARITY_CACHEDIR
chmod 777 $SINGULARITY_CACHEDIR

# Run the pipeline
echo "Starting LIANA + NicheNet pipeline"
echo "=================================="
echo "Input file:    $INPUT"
echo "Sender:        $SENDER"
echo "Receiver:      $RECEIVER"
echo "Output dir:    $OUTDIR"
echo "Profile:       $PROFILE"
echo "=================================="

# Run the pipeline with verbose logging enabled
NXF_VER=21.10.6 nextflow run main.nf \
  --input "$INPUT" \
  --sender_celltype "$SENDER" \
  --receiver_celltype "$RECEIVER" \
  --outdir "$OUTDIR" \
  -profile "$PROFILE" \
  -resume \
  -with-trace \
  -with-report \
  -ansi-log false

# Check exit status
if [[ $? -eq 0 ]]; then
  echo "Pipeline completed successfully!"
  echo "Results are available in: $OUTDIR"
else
  echo "Pipeline failed. Please check the logs for details."
  exit 1
fi
