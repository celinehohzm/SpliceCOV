#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ---------------------------
# Self-timing wrapper (GNU time)
# ---------------------------
_find_time_cmd() {
  if command -v /usr/bin/time >/dev/null 2>&1; then
    echo "/usr/bin/time"
  elif command -v gtime >/dev/null 2>&1; then
    echo "gtime"
  else
    echo ""
  fi
}

if [[ -z "${SPLICECOV_WRAPPED:-}" ]]; then
  _time_cmd="$(_find_time_cmd)"
  if [[ -n "$_time_cmd" ]]; then
    _timefile="$(mktemp -t splicecov.time.XXXXXX)"
    set +e
    SPLICECOV_WRAPPED=1 "$_time_cmd" -v -o "$_timefile" bash "$0" "$@"
    _rc=$?
    set -e
    if [[ -f "$_timefile" ]]; then
      _elapsed_str="$(grep -E 'Elapsed \(wall clock\) time' "$_timefile" | awk -F': ' '{print $2}')"
      _maxrss_kb="$(grep -E 'Maximum resident set size' "$_timefile" | awk -F': ' '{gsub(/ /,"",$2); print $2}')"
      to_seconds() {
        local IFS=':'; read -r -a p <<< "$1"
        local h=0 m=0 s=0
        if   ((${#p[@]}==3)); then h=${p[0]}; m=${p[1]}; s=${p[2]}
        elif ((${#p[@]}==2)); then m=${p[0]}; s=${p[1]}
        else                    s=${p[0]}; fi
        s=${s%.*}
        echo $((10#$h*3600 + 10#$m*60 + 10#$s))
      }
      _elapsed_sec="$(to_seconds "${_elapsed_str:-0}")"
      if [[ -n "${_maxrss_kb:-}" ]]; then
        _maxrss_mb=$(( (10#${_maxrss_kb} + 1023) / 1024 ))
        printf "[%(%F %T)T] Run summary: wall=%s (%ss), MaxRSS=%s KB (~%s MB)\n" -1 "${_elapsed_str:-N/A}" "${_elapsed_sec:-N/A}" "${_maxrss_kb}" "${_maxrss_mb}" >&2
      else
        printf "[%(%F %T)T] Run summary: wall=%s (%ss). (Install GNU time to report MaxRSS)\n" -1 "${_elapsed_str:-N/A}" "${_elapsed_sec:-N/A}" >&2
      fi
      rm -f "$_timefile"
    else
      printf "[%(%F %T)T] Run summary unavailable (timing file missing).\n" -1 >&2
    fi
    exit $_rc
  else
    _t0=$(date +%s)
    set +e
    SPLICECOV_WRAPPED=1 bash "$0" "$@"
    _rc=$?
    set -e
    _t1=$(date +%s)
    _dt=$((_t1 - _t0))
    printf "[%(%F %T)T] Run summary: wall=%ss. (Install GNU time for MaxRSS)\n" -1 "$_dt" >&2
    exit $_rc
  fi
fi

# ---------------------------
# Usage / CLI
# ---------------------------
usage() {
  cat <<'USAGE'
Usage: splicecov.sh -j <input_tiebrush_junc> -c <input_tiebrush_bigwig> [-a <annotation_gtf>] [-n <start_step>]

Required:
  -j <file> : input TieBrush junction file 
  -c <file> : input coverage BigWig file 

Optional:
  -a <file> : input annotation (GTF). If provided, annotation-dependent steps
              (building introns/unique splice sites and evaluation) will run.
  -n <int>  : step number to start from (default: 1).
              The pipeline has numbered steps (1 = sort junctions,
              2 = add coverage signal, 3 = build reference introns, …).
              By default, all steps run 1 → end. Use -n to resume later.

Notes:
  • All outputs are written to the "out/" directory.
  • If -a is omitted, evaluation steps are skipped automatically.

USAGE
  exit 1
}

# Default start step
start_step=1
input_annotation=""

# Parse arguments
while getopts ":j:c:a:n:h" opt; do
  case $opt in
    j) input_tiebrush_junc="$OPTARG" ;;
    c) input_tiebrush_bigwig="$OPTARG" ;;
    a) input_annotation="$OPTARG" ;;
    n) start_step="$OPTARG" ;;
    h) usage ;;
    \?) echo "Invalid option -$OPTARG" >&2; usage ;;
    :)  echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done

# Basic input validation
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

# ---------------------------
# Helpers: logging & checks
# ---------------------------
log() { printf "[%(%F %T)T] %s\n" -1 "$*" >&2; }
die() { echo "FATAL: $*" >&2; exit 1; }

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "'$1' not found in PATH"
}

need_file() {
  [[ -f "$1" ]] || die "Missing helper script: $1"
}

# ---------------------------
# Resolve script & helpers dir
# ---------------------------
this_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

helpers_dir="${SPLICECOV_HELPERS_DIR:-}"
if [[ -z "$helpers_dir" ]]; then
  for cand in \
      "$this_dir" \
      "$this_dir/scripts" \
      "$this_dir/../scripts"
  do
    if [[ -d "$cand" ]]; then
      helpers_dir="$cand"
      break
    fi
  done
fi
[[ -z "$helpers_dir" ]] && die "Could not locate helpers dir. Set SPLICECOV_HELPERS_DIR."

# List helper scripts we call
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
  "process_round1_splicecov_ptf_to_processed_juncs.py" # may be used with annotation workflows
)

# ---------------------------
# Preflight checks
# ---------------------------
log "Preflight: checking tools..."
need_cmd python
need_cmd awk
need_cmd sort
need_cmd bigWigToBedGraph

# Only needed when annotation is present
if $anno_present; then
  need_cmd comm
fi

# Verify helper scripts exist (conditional)
for s in "${common_helpers[@]}"; do
  need_file "${helpers_dir}/${s}"
done
if $anno_present; then
  for s in "${anno_helpers[@]}"; do
    need_file "${helpers_dir}/${s}"
  done
fi

# ---------------------------
# Outputs & paths
# ---------------------------
base_name="$(basename "$input_tiebrush_junc" .txt)"
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

# Converted BedGraph path (generated from the BigWig before Step 8)
converted_bedgraph="$outdir/${base_name}.bw.bedGraph"

# ---------------------------
# Pipeline
# ---------------------------
if [ "$start_step" -le 1 ]; then
  log "Step 1a: Sorting junctions by chr,start,end (header preserved)..."
  sorted_junc="$outdir/${base_name}.sorted.bed"

  sort_args=(-k1,1 -k2,2n -k3,3n)
  if LC_ALL=C sort --help 2>&1 | grep -q -- '--parallel'; then
    cpus="$(getconf _NPROCESSORS_ONLN 2>/dev/null || nproc 2>/dev/null || echo 1)"
    sort_args+=(--parallel="$cpus")
  fi
  if LC_ALL=C sort --help 2>&1 | grep -q -- '-T'; then
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
  python "${helpers_dir}/process_tiebrush_round1_juncs_splicecov.py" "$input_tiebrush_bigwig" "$processed_junc" > "$processed_junc_bundle"
fi

# ---------------------------
# Annotation-dependent block (Step 3 & Step 6)
# ---------------------------
if $anno_present && [ "$start_step" -le 3 ]; then
  log "Step 3: Building reference introns & unique splice sites (annotation provided)..."
  python "${helpers_dir}/gtf_to_intron_bed.py" "$input_annotation" "$annotation_introns"
  awk '{print $1, $2; print $1, $3}' "$annotation_introns" \
    | sort -k1,1 -k2,2n | uniq > "$annotation_uniquess"
  awk '{print $1, $2; print $1, $3+1}' "$input_tiebrush_junc" \
    | sort -k1,1 -k2,2n | uniq > "$tiebrush_uniquess"
  comm -12 <(sort "$annotation_uniquess") <(sort "$tiebrush_uniquess") \
    | sort -k1,1 -k2,2n | uniq > "$tiebrush_uniquess_in_annotation"
fi

if [ "$start_step" -le 4 ]; then
  log "Step 4: LightGBM scoring (junctions)..."
  python "${helpers_dir}/LightGBM_no_normscale.py" \
    -i "$processed_junc_bundle" -o "$processed_junc_bundle_w_scores"
fi

if [ "$start_step" -le 5 ]; then
  log "Step 5: Filtering score-positive junctions..."
  awk '$NF==1' "$processed_junc_bundle_w_scores" > "$processed_junc_bundle_w_scores_scpositive"
fi

# Step 6 depends on annotation
if $anno_present && [ "$start_step" -le 6 ]; then
  log "Step 6: Evaluation (junctions vs annotation splice sites)..."
  python "${helpers_dir}/evaluate_ptf_new.py" "$processed_junc_bundle_w_scores_scpositive" "$tiebrush_uniquess_in_annotation"
elif ! $anno_present && [ "$start_step" -le 6 ]; then
  log "Step 6: Skipped (no annotation provided)."
fi

if [ "$start_step" -le 7 ]; then
  log "Step 7: Emit PTF (junctions)..."
  awk 'BEGIN{OFS="\t"} { print $1, $2, $5, $12 }' \
    "$processed_junc_bundle_w_scores_scpositive" \
    > "$processed_junc_bundle_w_scores_scpositive_ptf"
fi

# ---- Round-2 TSSTES path ----
# Optional conversion to round1_processed_junc (kept commented as in original)
# if [ "$start_step" -le 8 ]; then
#   log "Step 7b: Convert PTF -> processed junctions (round 1)..."
#   python "${helpers_dir}/process_round1_splicecov_ptf_to_processed_juncs.py" "$processed_junc_bundle_w_scores_scpositive" "$processed_junc" "$round1_processed_junc"
# fi

# Step 8: Convert BigWig -> BedGraph then run round-2 re-processing
if [ "$start_step" -le 8 ]; then
  log "Step 8a: Converting BigWig -> BedGraph for round 2..."
  bigWigToBedGraph "$input_tiebrush_bigwig" "$converted_bedgraph"

  log "Step 8b: Re-processing original bedGraph for round 2..."
  # Note: round1_processed_junc may come from prior step (7b) if enabled; if not present,
  # the original script referenced it anyway. We proceed and let the helper decide or error.
  "${helpers_dir}/process_tiebrush_original.pl" \
    "$converted_bedgraph" "$round1_processed_junc" \
    > "$round2_processed_bundles"
fi

if [ "$start_step" -le 9 ]; then
  log "Step 9: Computing TSSTES metrics (round 2)..."
  python "${helpers_dir}/compute_round2_tsstes_metrics.py" \
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
  python "${helpers_dir}/LightGBM_tss.py" \
    -i "$round2_processed_bundles_w_metrics_tsstes_ptf" \
    -o "$round2_processed_bundles_w_metrics_tsstes_ptf_w_scores"
fi

if [ "$start_step" -le 13 ]; then
  log "Step 13: Filter TSSTES score-positive..."
  awk '$NF==1.0' \
    "$round2_processed_bundles_w_metrics_tsstes_ptf_w_scores" \
    > "$round2_processed_bundles_w_metrics_tsstes_ptf_w_scores_scpositive"
fi

# Step 14 depends on annotation
if $anno_present && [ "$start_step" -le 14 ]; then
  log "Step 14: Final TSSTES evaluation (annotation provided)..."
  python "${helpers_dir}/evaluate_TSSTES_new.py" \
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
