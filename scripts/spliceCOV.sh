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
Usage: splicecov.sh -j <input_tiebrush_junc> -c <input_tiebrush_bigwig> [-a <annotation_gtf>] [-b <basename>] [-s <threshold>]

Required:
  -j <file> : input TieBrush junction file 
  -c <file> : input coverage BigWig file 

Optional:
  -a <file> : input annotation (GTF). If provided, annotation-dependent steps
              (building introns/unique splice sites and evaluation) will run.
  -b <str>  : basename to use for ALL output files; overrides the default from -j.
  -s <num>  : probability threshold in [0,1] to pass to LightGBM scoring scripts (junctions & TSSTES).
              If omitted, those scripts use their own default (0.4).

Outputs (only these remain in out/):
  <basename>.jscore.txt
  <basename>.tsstes.scores.txt
  <basename>.combined.ptf
USAGE
  exit 1
}

# ---------------------------
# Defaults & arg parsing
# ---------------------------
input_annotation=""
basename_arg=""
score_arg=""

while getopts ":j:c:a:b:s:h" opt; do
  case $opt in
    j) input_tiebrush_junc="$OPTARG" ;;
    c) input_tiebrush_bigwig="$OPTARG" ;;
    a) input_annotation="$OPTARG" ;;
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

# >>> Models dir
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
need_cmd python3
need_cmd awk
need_cmd sort
need_cmd bigWigToBedGraph
$anno_present && need_cmd comm || true
for s in "${common_helpers[@]}"; do need_file "${helpers_dir}/${s}"; done
if $anno_present; then for s in "${anno_helpers[@]}"; do need_file "${helpers_dir}/${s}"]; done; fi

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

# Final deliverables (only these remain)
jscore_out="$outdir/${base_name}.jscore.txt"
tsstes_scores_out="$outdir/${base_name}.tsstes.scores.txt"
combined_out="$outdir/${base_name}.combined.ptf"

# LightGBM threshold flags
declare -a score_flags=()
if [[ -n "$score_arg" ]]; then
  score_flags=(-s "$score_arg")
fi

# Temp workspace for intermediates
workdir="$(mktemp -d "${outdir}/.work.${base_name}.XXXX")"
cleanup() {
  if [[ -n "${KEEP_TEMP:-}" ]]; then
    log "KEEP_TEMP set; leaving workspace: $workdir"
  else
    rm -rf "$workdir"
  fi
}
trap cleanup EXIT

# Intermediates -> workdir
sorted_junc="$workdir/${base_name}.sorted.bed"
processed_junc="$workdir/${base_name}.jproc.txt"
processed_junc_bundle="$workdir/${base_name}.jbund.txt"
jpos_source="$workdir/${base_name}.jpos.txt"          # score-positive junctions (pre-PTF)
jpos_ptf_tmp="$workdir/${base_name}.jpos.ptf"         # temp PTF for junctions (not persisted)
tsstes_pos_tmp="$workdir/${base_name}.tsstes.pos.txt" # temp positives (not persisted)

# Annotation-dependent intermediates
annotation_introns="$workdir/ann.introns.bed"
annotation_uniquess="$workdir/${base_name}.ann.uss.bed"
tiebrush_uniquess="$workdir/${base_name}.base.uss.bed"
tiebrush_uniquess_in_annotation="$workdir/${base_name}.base.uss.in_ann.bed"
tiebrush_unique_tsstes_in_annotation="$workdir/${base_name}.base.tsstes.in_ann.bed" # (may remain unused)

# Round-2 / TSSTES intermediates
round1_processed_junc="$workdir/${base_name}.r1.jproc.txt"
round2_processed_bundles="$workdir/${base_name}.bund.txt"
round2_processed_bundles_w_metrics="$workdir/${base_name}.r2.metrics.txt"
round2_processed_bundles_w_metrics_ptf="$workdir/${base_name}.r2.metrics.ptf"
round2_processed_bundles_w_metrics_tsstes_ptf="$workdir/${base_name}.tsstes.ptf"
converted_bedgraph="$workdir/${base_name}.bw.bedGraph"

# ---------------------------
# Pipeline (always from start)
# ---------------------------
log "Step 1a: Sorting junctions by chr,start,end (header preserved)..."
sort_args=(-k1,1 -k2,2n -k3,3n)
if LC_ALL=C sort --help 2>/dev/null | grep -q -- '--parallel'; then
  cpus="$(getconf _NPROCESSORS_ONLN 2>/dev/null || nproc 2>/dev/null || echo 1)"
  sort_args+=(--parallel="$cpus")
fi
if LC_ALL=C sort --help 2>/dev/null | grep -q -- '-T'; then
  sort_args+=(-T "${TMPDIR:-$workdir}")
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

log "Step 2: Adding bigWig signal..."
python3 "${helpers_dir}/process_tiebrush_round1_juncs_splicecov.py" \
  "$input_tiebrush_bigwig" "$processed_junc" > "$processed_junc_bundle"

# ---------------------------
# Annotation-dependent block — results not persisted
# ---------------------------
if $anno_present; then
  log "Step 3: Building reference introns & unique splice sites (annotation provided)..."
  python3 "${helpers_dir}/gtf_to_intron_bed.py" "$input_annotation" "$annotation_introns"
  awk '{print $1, $2; print $1, $3}' "$annotation_introns" \
    | sort -k1,1 -k2,2n | uniq > "$annotation_uniquess"
  awk '{print $1, $2; print $1, $3+1}' "$input_tiebrush_junc" \
    | sort -k1,1 -k2,2n | uniq > "$tiebrush_uniquess"
  comm -12 <(sort "$annotation_uniquess") <(sort "$tiebrush_uniquess") \
    | sort -k1,1 -k2,2n | uniq > "$tiebrush_uniquess_in_annotation"
fi

log "Step 4: LightGBM scoring (junctions) -> ${jscore_out}"
python3 "${helpers_dir}/LightGBM_no_normscale.py" \
  -i "$processed_junc_bundle" \
  -o "$jscore_out" \
  ${score_flags[@]+"${score_flags[@]}"}

log "Step 5: Filtering score-positive junctions..."
awk '$NF==1 || $NF==1.0' "$jscore_out" > "$jpos_source"

if $anno_present; then
  log "Step 6: (Optional) Evaluation vs annotation (not persisted)..."
  if [[ -s "$jpos_source" && -s "$tiebrush_uniquess_in_annotation" ]]; then
    set +e
    python "${helpers_dir}/evaluate_ptf_new.py" \
      "$jpos_source" \
      "$tiebrush_uniquess_in_annotation" || true
    set -e
  else
    log "Step 6: Skipped (no positives or no annotation splice sites)."
  fi
fi

log "Step 7: Emit PTF (junctions) -> temp (not persisted)"
awk 'BEGIN{OFS="\t"} { print $1, $2, $5, $12 }' \
  "$jpos_source" > "$jpos_ptf_tmp"

log "Step 8a: Converting BigWig -> BedGraph for round 2..."
bigWigToBedGraph "$input_tiebrush_bigwig" "$converted_bedgraph"

log "Step 8b: Re-processing original bedGraph for round 2..."
"${helpers_dir}/process_tiebrush_original.pl" \
  "$converted_bedgraph" "$round1_processed_junc" \
  > "$round2_processed_bundles"

log "Step 9: Computing TSSTES metrics (round 2)..."
python3 "${helpers_dir}/compute_round2_tsstes_metrics.py" \
  "$round2_processed_bundles" "$input_tiebrush_bigwig" \
  > "$round2_processed_bundles_w_metrics"

log "Step 10: Bundles -> PTF (round 2 TSSTES)..."
"${helpers_dir}/splicecov_bundle2ptf.pl" \
  "$round2_processed_bundles_w_metrics" \
  > "$round2_processed_bundles_w_metrics_ptf"

log "Step 11: Extract TSS/CPAS..."
awk '($4=="TSS" || $4=="CPAS")' \
  "$round2_processed_bundles_w_metrics_ptf" \
  > "$round2_processed_bundles_w_metrics_tsstes_ptf"

log "Step 12: LightGBM scoring (TSSTES) -> ${tsstes_scores_out}"
python3 "${helpers_dir}/LightGBM_tss.py" \
  -i "$round2_processed_bundles_w_metrics_tsstes_ptf" \
  -o "$tsstes_scores_out" \
  ${score_flags[@]+"${score_flags[@]}"}

log "Step 13: Filter TSSTES score-positive -> temp (not persisted)"
awk '$NF==1 || $NF==1.0' \
  "$tsstes_scores_out" > "$tsstes_pos_tmp"

if $anno_present; then
  log "Step 14: (Optional) TSSTES evaluation vs annotation (not persisted)..."
  if [[ -s "$tsstes_pos_tmp" && -s "$tiebrush_unique_tsstes_in_annotation" ]]; then
    set +e
    python "${helpers_dir}/evaluate_TSSTES_new.py" \
      "$tsstes_pos_tmp" \
      "$tiebrush_unique_tsstes_in_annotation" \
      > /dev/null || true
    set -e
  else
    log "Step 14: Skipped (no TSSTES positives or no annotation TSSTES reference)."
  fi
fi

log "Step 15: Combine ptfs -> ${combined_out}"
"${helpers_dir}/combine_ptfs.sh" \
  "$jpos_source" \
  "$tsstes_pos_tmp" \
  > "$combined_out"

log "Final outputs:"
ls -lh "$jscore_out" "$tsstes_scores_out" "$combined_out" || true

log "Cleaning up intermediates..."
# handled by trap cleanup on $workdir

log "spliceCOV complete!"
