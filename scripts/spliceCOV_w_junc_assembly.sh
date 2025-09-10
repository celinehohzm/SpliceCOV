#!/bin/bash

set -e  # Exit on any error

# Function to display usage
usage() {
  echo "Usage: $0 -j <input_tiebrush_junc> -c <input_tiebrush_bigwig> -b <input_tiebrush_bedgraph> -a <input_annotation>"
  exit 1
}

# Parse arguments
while getopts ":j:c:b:a:" opt; do
  case $opt in
    j) input_tiebrush_junc="$OPTARG" ;;
    c) input_tiebrush_bigwig="$OPTARG" ;;
    b) input_tiebrush_bedgraph="$OPTARG" ;;
    a) input_annotation="$OPTARG" ;;
    \?) echo "Invalid option -$OPTARG" >&2; usage ;;
    :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done

# Check that all inputs are provided
if [ -z "$input_tiebrush_junc" ] || [ -z "$input_tiebrush_bigwig" ] || [ -z "$input_tiebrush_bedgraph" ] || [ -z "$input_annotation" ]; then
  usage
fi

# Generate output filenames based on input junction file name
base_name=$(basename "$input_tiebrush_junc" .txt)
outdir="out"
mkdir -p "$outdir"

processed_junc="$outdir/${base_name}_processed.txt"
processed_junc_bundle="$outdir/${base_name}_bundle.txt"
processed_junc_bundle_w_scores="$outdir/${base_name}_bundle_scores.txt"
processed_junc_bundle_w_scores_scpositive="$outdir/${base_name}_bundle_scores_scpositive.txt"
processed_junc_bundle_w_scores_scpositive_ptf="$outdir/${base_name}_bundle_scores_scpositive.ptf"

annotation_introns="$outdir/annotation_introns.bed"
annotation_uniquess="$outdir/${base_name}_annotation_uniquess.bed"
tiebrush_uniquess="$outdir/${base_name}_tiebrush_uniquess.bed"
tiebrush_uniquess_in_annotation="$outdir/${base_name}_tiebrush_uniquess_in_annotation.bed"

# Variables for round-2 steps
round1_processed_junc="$outdir/${base_name}_round1_processed_juncs.txt"

# Directory containing the scripts
script_dir="/ccb/salz1/choh1/spliceCov/scripts"

# Step 1
echo "Step 1: Processing junctions..."
"${script_dir}/process_junctions_perc.pl" "$input_tiebrush_junc" > "$processed_junc"

# Step 2
echo "Step 2: Adding bigwig signal..."
python "${script_dir}/process_tiebrush_round1_juncs_splicecov.py" \
  "$input_tiebrush_bigwig" "$processed_junc" > "$processed_junc_bundle"

# Step 3
echo "Step 3: Running LightGBM scoring..."
python "${script_dir}/LightGBM.py" -i "$processed_junc_bundle" -o "$processed_junc_bundle_w_scores"

# Step 4
echo "Step 4: Creating reference file for evaluation (tiebrush in $input_annotation)..."
python "${script_dir}/gtf_to_intron_bed.py" "$input_annotation" "$annotation_introns"
awk '{print $1, $2; print $1, $3}' "$annotation_introns" \
  | sort -k1,1 -k2,2n | uniq > "$annotation_uniquess"
awk '{print $1, $2; print $1, $3+1}' "$input_tiebrush_junc" \
  | sort -k1,1 -k2,2n | uniq > "$tiebrush_uniquess"
comm -12 <(sort "$annotation_uniquess") <(sort "$tiebrush_uniquess") \
  | sort -k1,1 -k2,2n | uniq > "$tiebrush_uniquess_in_annotation"

# Step 5
echo "Step 5: Filtering score-positive junctions..."
awk '$NF==1' "$processed_junc_bundle_w_scores" > "$processed_junc_bundle_w_scores_scpositive"

# Step 6
echo "Step 6: Final evaluation..."
python "${script_dir}/evaluate_ptf_new.py" \
  "$processed_junc_bundle_w_scores_scpositive" "$tiebrush_uniquess_in_annotation"

# Step 7
echo "Step 7: Gets ptf file"
awk 'BEGIN{OFS="\t"} { print $1, $2, $5, $12 }' "$processed_junc_bundle_w_scores_scpositive" > "$processed_junc_bundle_w_scores_scpositive_ptf"


echo "spliceCOV complete!"
