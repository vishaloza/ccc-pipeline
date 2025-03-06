#!/bin/bash
# Helper script to run the LIANA + NicheNet pipeline

# Default values
INPUT=""
SENDER=""
RECEIVER=""
OUTDIR="results"
PROFILE="conda"

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
  echo "Usage: ./run.sh --input <h5ad_file> --sender <sender_celltype> --receiver <receiver_celltype> [--outdir <output_dir>] [--profile <nextflow_profile>]"
  exit 1
fi

# Check if file exists
if [[ ! -f "$INPUT" ]]; then
  echo "Error: Input file not found: $INPUT"
  exit 1
fi

# Run Nextflow pipeline
echo "Starting LIANA + NicheNet pipeline"
echo "=================================="
echo "Input file:    $INPUT"
echo "Sender:        $SENDER"
echo "Receiver:      $RECEIVER"
echo "Output dir:    $OUTDIR"
echo "Profile:       $PROFILE"
echo "=================================="

# Create output directory
mkdir -p "$OUTDIR"

# Run the pipeline
nextflow run main.nf \
  --input "$INPUT" \
  --sender_celltype "$SENDER" \
  --receiver_celltype "$RECEIVER" \
  --outdir "$OUTDIR" \
  -profile "$PROFILE" \
  -resume

# Check exit status
if [[ $? -eq 0 ]]; then
  echo "Pipeline completed successfully!"
  echo "Results are available in: $OUTDIR"
else
  echo "Pipeline failed. Please check the logs for details."
  exit 1
fi
