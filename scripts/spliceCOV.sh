#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ---------------------------
# Portable logging & helpers
# ---------------------------
ts() { date '+%F %T'; }
log() { printf '[%s] %s\n' "$(ts)" "$*" >&2; }
die() { echo "FATAL: $*" >&2; exit 1; }

need_cmd()  { command -v "$1" >/dev/null 2>&1 || die "'$1' not found in PATH"; }
need_file() { [[ -f "$1" ]] || die "Missing helper script: $1"; }

# ---------------------------
# Usage / CLI
# ---------------------------
usage() {
  cat <<'USAGE'
Usage: splicecov.sh -j <input_tiebrush_junc> -c <input_tiebrush_bigwig> [-a <annotation_gtf>] [-n <start_step>] [-b <basename>] [-s <threshold>]

Required:
  -j <file> : input TieBrush junction file 
  -c <file> : input coverage BigWig file 

Optional:
  -a <file> : input annotation (GTF). If provided, annotation-dependent steps
              (building introns/unique splice sites and evaluation) will run.
  -n <int>  : step number to start from (default: 1). 1 = sort junctions, 2 = add coverage, ...
  -b <str>  : basename to use for ALL output files; overrides the default from -j.
  -s <num>  : probability threshold in [0,1] to pass to LightGBM scoring scripts (steps 4 & 12).
              If omitted, those scripts use their own default (0.4).

Notes:
  • All outputs are written to the "out/" directory.
  • If -a is omitted, evaluation steps are skipped automatically.
USAGE
  exit 1
}

# ---------------------------
# Defaults & arg parsing
# ---------------------------
start_step=1
input_annotation=""
basename_arg=""
score_arg=""

while getopts ":j:c:a:n:b:s:h" opt; do
  case $opt in
    j) input_tiebrush_junc="$OPTARG" ;;
    c) input_tiebrush_bigwig="$OPTARG" ;;
    a) input_annotation="$OPTARG" ;;
    n) start_step="$OPTARG" ;;
    b) basename_arg="$OPTARG" ;;
    s) score_arg="$OPTARG" ;;
    h) usage ;;
    \?) echo "Invalid option -$OPTARG" >&2; usage ;;
    :)  echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done

# ---------------------------
# Validate inputs
# ---------------------------
[[ -z "${input_tiebrush_junc:-}" || -z "${input_tiebrush_bigwig:-}" ]] && usage
[[ ! -f "$input_tiebrush_junc" ]]   && { echo "ERROR: Junction file not found: $input_tiebrush_junc" >&2; exit 2; }
[[ ! -f "$input_tiebrush_bigwig" ]] && { echo "ERROR: BigWig file not found: $input_tiebrush_bigwig" >&2; exit 2; }

anno_present=false
if [[ -n "$input_annotation" ]]; then
  [[ ! -f "$input_annotation" ]] && { echo "ERROR: Annotation GTF not found: $input_annotation" >&2; exit 2; }
  anno_present=true
fi

if ! [[ "$start_step" =~ ^[0-9]+$ ]] || [ "$start_step" -lt 1 ]; then
  echo "ERROR: -n requires a positive integer." >&2; exit 2
fi

if [[ -n "$basename_arg" ]]; then
  if [[ "$basename_arg" =~ [/\ ] ]]; then
    echo "ERROR: -b basename must not contain slashes or spaces." >&2; exit 2
  fi
fi

# Validate -s if provided: allow 0, 1, or decimals between
if [[ -n "$score_arg" ]]; then
  if ! [[ "$score_arg" =~ ^(0(\.[0-9]+)?|1(\.0+)?)$ ]]; then
    echo "ERROR: -s must be a number in [0,1], got '$score_arg'." >&2
    exit 2
  fi
fi

# ---------------------------
# Resolve script & helpers dir
# ---------------------------
this_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

helpers_dir="${SPLICECOV_HELPERS_DIR:-}"
if [[ -z "$helpers_dir" ]]; then
  for cand in "$this_dir" "$this_dir/scripts" "$this_dir/../scripts"; do
    if [[ -d "$cand" ]]; then
      helpers_dir="$cand"; break
    fi
  done
fi
[[ -z "$helpers_dir" ]] && die "Could not locate helpers dir. Set SPLICECOV_HELPERS_DIR."

# >>> NEW: resolve model directory once (env override supported)
MODEL_DIR="${SPLICECOV_MODEL_DIR:-${helpers_dir%/}/model_output}"
export SPLICECOV_MODEL_DIR="$MODEL_DIR"
log "Models dir: $MODEL_DIR"
# <<<

# Helper scripts we call
common_helpers=(
  "process_junctions_perc.pl"
  "process_tiebrush_round1_juncs_splicecov.py"
  "LightGBM_no_normscale.py"
  "process_tiebrush_original.pl"
  "compute_round2_tsstes_metrics.py"
  "splicecov_bundle2ptf.pl"
  "LightGBM_tss.py"
  "combine_ptfs.sh"
)
anno_helpers=(
  "gtf_to_intron_bed.py"
  "evaluate_ptf_new.py"
  "evaluate_TSSTES_new.py"
  "process_round1_splicecov_ptf_to_processed_juncs.py"
)

# ---------------------------
# Preflight checks
# ---------------------------
log "Preflight: checking tools..."
need_cmd python3             # changed to python3
need_cmd awk
need_cmd sort
need_cmd bigWigToBedGraph
$anno_present && need_cmd comm || true
for s in "${common_helpers[@]}"; do need_file "${helpers_dir}/${s}"; done
if $anno_present; then for s in "${anno_helpers[@]}"; do need_file "${helpers_dir}/${s}"; done; fi

# ---------------------------
# Outputs & paths
# ---------------------------
if [[ -n "$basename_arg" ]]; then
  base_name="$basename_arg"
else
  f="${input_tiebrush_junc##*/}"
  if [[ "$f" == *.* ]]; then base_name="${f%.*}"; else base_name="$f"; fi
fi

outdir="out"
mkdir -p "$outdir"

processed_junc="$outdir/${base_name}.jproc.txt"
processed_junc_bundle="$outdir/${base_name}.jbund.txt"
processed_junc_bundle_w_scores="$outdir/${base_name}.jscore.txt"
processed_junc_bundle_w_scores_scpositive="$outdir/${base_name}.jpos.txt"
processed_junc_bundle_w_scores_scpositive_ptf="$outdir/${base_name}.jpos.ptf"

# Annotation-dependent outputs
annotation_introns="$outdir/ann.introns.bed"
annotation_uniquess="$outdir/${base_name}.ann.uss.bed"
tiebrush_uniquess="$outdir/${base_name}.base.uss.bed"
tiebrush_uniquess_in_annotation="$outdir/${base_name}.base.uss.in_ann.bed"
tiebrush_unique_tsstes_in_annotation="$outdir/${base_name}.base.tsstes.in_ann.bed"

# Round-2 / TSSTES outputs
round1_processed_junc="$outdir/${base_name}.r1.jproc.txt"
round2_processed_bundles="$outdir/${base_name}.r2.bund.txt"
round2_processed_bundles_w_metrics="$outdir/${base_name}.r2.metrics.txt"
round2_processed_bundles_w_metrics_ptf="$outdir/${base_name}.r2.metrics.ptf"
round2_processed_bundles_w_metrics_tsstes_ptf="$outdir/${base_name}r2.tsstes.ptf"
round2_processed_bundles_w_metrics_tsstes_ptf_w_scores="$outdir/${base_name}.r2.tsstes.scores.txt"
round2_processed_bundles_w_metrics_tsstes_ptf_w_scores_scpositive="$outdir/${base_name}.r2.tsstes.pos.txt"
round2_processed_bundles_w_metrics_tsstes_ptf_w_scores_scpositive_eval="$outdir/${base_name}.r2.tsstes.pos.eval.txt"

converted_bedgraph="$outdir/${base_name}.bw.bedGraph"

# LightGBM threshold flags
declare -a score_flags=()
if [[ -n "$score_arg" ]]; then
  score_flags=(-s "$score_arg")
fi

# ---------------------------
# Pipeline
# ---------------------------
if [ "$start_step" -le 1 ]; then
  log "Step 1a: Sorting junctions by chr,start,end (header preserved)..."
  sorted_junc="$outdir/${base_name}.sorted.bed"

  sort_args=(-k1,1 -k2,2n -k3,3n)
  if LC_ALL=C sort --help 2>/dev/null | grep -q -- '--parallel'; then
    cpus="$(getconf _NPROCESSORS_ONLN 2>/dev/null || nproc 2>/dev/null || echo 1)"
    sort_args+=(--parallel="$cpus")
  fi
  if LC_ALL=C sort --help 2>/dev/null | grep -q -- '-T'; then
    sort_args+=(-T "${TMPDIR:-$outdir}")
  fi

  first_line="$(head -n1 "$input_tiebrush_junc")"
  {
    if [[ "$first_line" =~ ^(track|#) ]]; then
      printf '%s\n' "$first_line"
      tail -n +2 "$input_tiebrush_junc" | LC_ALL=C sort "${sort_args[@]}"
    else
      printf 'track name=junctions\n'
      LC_ALL=C sort "${sort_args[@]}" "$input_tiebrush_junc"
    fi
  } > "$sorted_junc"

  log "Step 1b: Processing junctions (sorted input)..."
  "${helpers_dir}/process_junctions_perc.pl" "$sorted_junc" > "$processed_junc"
fi

if [ "$start_step" -le 2 ]; then
  log "Step 2: Adding bigWig signal..."
  python3 "${helpers_dir}/process_tiebrush_round1_juncs_splicecov.py" \
    "$input_tiebrush_bigwig" "$processed_junc" > "$processed_junc_bundle"
fi

# ---------------------------
# Annotation-dependent block (Step 3 & Step 6)
# ---------------------------
if $anno_present && [ "$start_step" -le 3 ]; then
  log "Step 3: Building reference introns & unique splice sites (annotation provided)..."
  python3 "${helpers_dir}/gtf_to_intron_bed.py" "$input_annotation" "$annotation_introns"
  awk '{print $1, $2; print $1, $3}' "$annotation_introns" \
    | sort -k1,1 -k2,2n | uniq > "$annotation_uniquess"
  awk '{print $1, $2; print $1, $3+1}' "$input_tiebrush_junc" \
    | sort -k1,1 -k2,2n | uniq > "$tiebrush_uniquess"
  comm -12 <(sort "$annotation_uniquess") <(sort "$tiebrush_uniquess") \
    | sort -k1,1 -k2,2n | uniq > "$tiebrush_uniquess_in_annotation"
fi

if [ "$start_step" -le 4 ]; then
  log "Step 4: LightGBM scoring (junctions)..."
  python3 "${helpers_dir}/LightGBM_no_normscale.py" \
    -i "$processed_junc_bundle" \
    -o "$processed_junc_bundle_w_scores" \
    ${score_flags[@]+"${score_flags[@]}"}
fi

if [ "$start_step" -le 5 ]; then
  log "Step 5: Filtering score-positive junctions..."
  awk '$NF==1' "$processed_junc_bundle_w_scores" > "$processed_junc_bundle_w_scores_scpositive"
fi

# Step 6 depends on annotation
if $anno_present && [ "$start_step" -le 6 ]; then
  log "Step 6: Evaluation (junctions vs annotation splice sites)..."
  python3 "${helpers_dir}/evaluate_ptf_new.py" \
    "$processed_junc_bundle_w_scores_scpositive" "$tiebrush_uniquess_in_annotation"
elif ! $anno_present && [ "$start_step" -le 6 ]; then
  log "Step 6: Skipped (no annotation provided)."
fi

if [ "$start_step" -le 7 ]; then
  log "Step 7: Emit PTF (junctions)..."
  awk 'BEGIN{OFS="\t"} { print $1, $2, $5, $12 }' \
    "$processed_junc_bundle_w_scores_scpositive" \
    > "$processed_junc_bundle_w_scores_scpositive_ptf"
fi

if [ "$start_step" -le 8 ]; then
  log "Step 8a: Converting BigWig -> BedGraph for round 2..."
  bigWigToBedGraph "$input_tiebrush_bigwig" "$converted_bedgraph"

  log "Step 8b: Re-processing original bedGraph for round 2..."
  "${helpers_dir}/process_tiebrush_original.pl" \
    "$converted_bedgraph" "$round1_processed_junc" \
    > "$round2_processed_bundles"
fi

if [ "$start_step" -le 9 ]; then
  log "Step 9: Computing TSSTES metrics (round 2)..."
  python3 "${helpers_dir}/compute_round2_tsstes_metrics.py" \
    "$round2_processed_bundles" "$input_tiebrush_bigwig" \
    > "$round2_processed_bundles_w_metrics"
fi

if [ "$start_step" -le 10 ]; then
  log "Step 10: Bundles -> PTF (round 2 TSSTES)..."
  "${helpers_dir}/splicecov_bundle2ptf.pl" \
    "$round2_processed_bundles_w_metrics" \
    > "$round2_processed_bundles_w_metrics_ptf"
fi

if [ "$start_step" -le 11 ]; then
  log "Step 11: Extract TSS/CPAS..."
  awk '($4=="TSS" || $4=="CPAS")' \
    "$round2_processed_bundles_w_metrics_ptf" \
    > "$round2_processed_bundles_w_metrics_tsstes_ptf"
fi

if [ "$start_step" -le 12 ]; then
  log "Step 12: LightGBM scoring (TSSTES)..."
  python3 "${helpers_dir}/LightGBM_tss.py" \
    -i "$round2_processed_bundles_w_metrics_tsstes_ptf" \
    -o "$round2_processed_bundles_w_metrics_tsstes_ptf_w_scores" \
    ${score_flags[@]+"${score_flags[@]}"}
fi

if [ "$start_step" -le 13 ]; then
  log "Step 13: Filter TSSTES score-positive..."
  awk '$NF==1.0' \
    "$round2_processed_bundles_w_metrics_tsstes_ptf_w_scores" \
    > "$round2_processed_bundles_w_metrics_tsstes_ptf_w_scores_scpositive"
fi

if $anno_present && [ "$start_step" -le 14 ]; then
  log "Step 14: Final TSSTES evaluation (annotation provided)..."
  python3 "${helpers_dir}/evaluate_TSSTES_new.py" \
    "$round2_processed_bundles_w_metrics_tsstes_ptf_w_scores_scpositive" \
    "$tiebrush_unique_tsstes_in_annotation" \
    > "$round2_processed_bundles_w_metrics_tsstes_ptf_w_scores_scpositive_eval"
elif ! $anno_present && [ "$start_step" -le 14 ]; then
  log "Step 14: Skipped (no annotation provided)."
fi

if [ "$start_step" -le 15 ]; then
  log "Step 15: Combine ptfs..."
  "${helpers_dir}/combine_ptfs.sh" \
    "$processed_junc_bundle_w_scores_scpositive" \
    "$round2_processed_bundles_w_metrics_tsstes_ptf_w_scores_scpositive" \
    > "$outdir/${base_name}_combined.ptf"
fi

log "spliceCOV complete!"
