#!/usr/bin/env bash
# combine_sites.sh
# Usage: combine_sites.sh <bundle_scores_scpositive.txt> <tsstes_ptf_w_scores_scpositive.txt> [output.tsv]
# Output columns: chrom  pos  field3  field4
#   - from file1 (bundle scores): $1,$2,$5,$12   (chrom, pos, strand, JSTART/JEND)
#   - from file2 (TSSTES PTF):    $1,$2,$3,$4    (chrom, pos, ".", TSS/CPAS)
# Deduplicates identical rows.

set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 <bundle_scores_scpositive.txt> <tsstes_ptf_w_scores_scpositive.txt> [output.tsv]" >&2
  exit 1
fi

file1="$1"
file2="$2"
out="${3:-/dev/stdout}"

[[ -f "$file1" ]] || { echo "Missing: $file1" >&2; exit 2; }
[[ -f "$file2" ]] || { echo "Missing: $file2" >&2; exit 2; }

# Optional: tune sort tmp and threads
sort_tmp="${TMPDIR:-.}"
threads="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)"

{
  # First file: chrom, pos, strand, flag
  awk 'BEGIN{FS=OFS="\t"} !/^#/ && NF>=12 {print $1,$2,$5,$12}' "$file1"

  # Second file: chrom, pos, ".", TSS/CPAS
  awk 'BEGIN{FS=OFS="\t"} !/^#/ && NF>=4  {print $1,$2,$3,$4}'  "$file2"
} |
LC_ALL=C sort -T "$sort_tmp" --parallel="$threads" \
  -k1,1 -k2,2n -k3,3 -k4,4 \
  -u > "$out"
