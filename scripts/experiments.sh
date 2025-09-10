#!/usr/bin/env bash

#
# run_all_tissues.sh
#
# Usage:
#   ./run_all_tissues.sh
#
# This script iterates over a predefined list of tissues, constructs the “Capitalized”
# version of each (e.g. small_intestine → Small_Intestine), and then invokes spliceCOV.sh
# with the appropriate paths for each tissue.
#

set -euo pipefail

# List of tissues (lowercase, underscores)
tissues=(
  bladder
  small_intestine
  breast
  cervix_uteri
  colon
  fallopian_tube
  heart
  lung
  muscle
  nerve
  ovary
  pancreas
  pituitary
  prostate
  salivary_gland
  spleen
  stomach
  thyroid
  uterus
  vagina
)

# Base directories and script path
SPLICE_COV_SCRIPT="/ccb/salz1/choh1/spliceCov/scripts/spliceCOV.sh"
TIEBRUSH_RESULTS_BASE="/ccb/salz4-2/mpertea/GTEX_stringtie/chess3/tiebrush_results"
BIGWIG_BASE="/ccb/salz1/choh1/gtex_tiebrush_bigwig"
ANNOTATION_GTF="/home/choh1/spliceCOV/data/chess3.0.1.gtf"

# Loop over each tissue
for tissue in "${tissues[@]}"; do
  # Split on underscore to capitalize each part
  IFS='_' read -r -a parts <<< "$tissue"
  capitalized_parts=()
  for part in "${parts[@]}"; do
    first_char="${part:0:1}"
    rest_chars="${part:1}"
    capitalized_parts+=( "$(tr '[:lower:]' '[:upper:]' <<< "$first_char")${rest_chars}" )
  done
  capitalized_tissue="$(IFS=_; echo "${capitalized_parts[*]}")"

  # Construct paths
  junc_file="$TIEBRUSH_RESULTS_BASE/$capitalized_tissue/${capitalized_tissue}.def.junctions.bed"
  bigwig_file="$BIGWIG_BASE/$capitalized_tissue/${capitalized_tissue}.def.coverage.bigwig"
  bedgraph_file="$TIEBRUSH_RESULTS_BASE/$capitalized_tissue/${capitalized_tissue}.def.coverage.bedgraph"

  echo "Running spliceCOV for tissue: $capitalized_tissue"
  "$SPLICE_COV_SCRIPT" \
    -j "$junc_file" \
    -c "$bigwig_file" \
    -b "$bedgraph_file" \
    -a "$ANNOTATION_GTF"

  echo "Completed: $capitalized_tissue"
  echo "----------------------------------------"
done

echo "All tissues processed."
