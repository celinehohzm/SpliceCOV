#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

ts() { date '+%F %T'; }
log() { printf '[%s] %s\n' "$(ts)" "$*" >&2; }
die() { echo "FATAL: $*" >&2; exit 1; }

need_cmd()  { command -v "$1" >/dev/null 2>&1 || die "'$1' not found in PATH"; }
need_file() { [[ -f "$1" ]] || die "Missing helper script: $1"; }

usage() {
  cat <<'USAGE'
Usage:
  Full run:
    splicecov -j <input_tiebrush_junc> -c <input_tiebrush_bigwig> [-a <annotation_gtf>] [-b <basename>] [-s <threshold>]

  Eval-only (no pipeline; reads outputs from out/):
    splicecov -b <basename> -a <annotation_gtf>

Required for full run:
  -j <file> : input TieBrush junction file
  -c <file> : input coverage BigWig file

Required for eval-only:
  -b <str>  : basename used in out/<basename>.{jscore,tsstes.scores}.txt
  -a <file> : annotation GTF

Optional:
  -a <file> : annotation (GTF). If provided, evaluation will run *based on generated out/ files*.
  -b <str>  : basename to use for ALL output files; overrides default from -j.
  -s <num>  : LightGBM scoring threshold [0,1], default is 0.4.

Core outputs (written to out/ when running full pipeline):
  <basename>.jscore.txt
  <basename>.tsstes.scores.txt
  <basename>.combined.ptf

Evaluation outputs (written to out/ when -a is provided):
  <basename>.eval.junctions.txt
  <basename>.eval.tsstes.txt
USAGE
  exit 1
}

input_tiebrush_junc=""
input_tiebrush_bigwig=""
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

need_cmd python3
need_cmd awk
need_cmd sort

# Validate -s if provided
if [[ -n "$score_arg" ]]; then
  if ! [[ "$score_arg" =~ ^(0(\.[0-9]+)?|1(\.0+)?)$ ]]; then
    echo "ERROR: -s must be a number in [0,1], got '$score_arg'." >&2
    exit 2
  fi
fi

# Determine mode:
# - Full run if (-j and -c) provided
# - Eval-only if (-b and -a) provided and (-j/-c) not provided
full_run=false
eval_only=false

if [[ -n "$input_tiebrush_junc" || -n "$input_tiebrush_bigwig" ]]; then
  # If either is set, require both.
  [[ -z "$input_tiebrush_junc" || -z "$input_tiebrush_bigwig" ]] && usage
  full_run=true
else
  # No -j/-c means eval-only must be triggered via -b -a
  if [[ -n "$basename_arg" && -n "$input_annotation" ]]; then
    eval_only=true
  else
    usage
  fi
fi

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

# Resolve script & helpers dir
this_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

helpers_dir="${SPLICECOV_HELPERS_DIR:-}"
if [[ -z "$helpers_dir" ]]; then
  for cand in "$this_dir" "$this_dir/scripts" "$this_dir/../scripts"; do
    if [[ -d "$cand" ]]; then helpers_dir="$cand"; break; fi
  done
fi
[[ -z "$helpers_dir" ]] && die "Could not locate helpers dir. Set SPLICECOV_HELPERS_DIR."

MODEL_DIR="${SPLICECOV_MODEL_DIR:-${helpers_dir%/}/model_output}"
export SPLICECOV_MODEL_DIR="$MODEL_DIR"
log "Models dir: $MODEL_DIR"

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
)

log "Preflight: checking helper scripts..."
for s in "${common_helpers[@]}"; do need_file "${helpers_dir}/${s}"; done
if $anno_present; then for s in "${anno_helpers[@]}"; do need_file "${helpers_dir}/${s}"; done; fi

# Output paths
outdir="out"
mkdir -p "$outdir"

# Basename resolution
if [[ -n "$basename_arg" ]]; then
  base_name="$basename_arg"
else
  f="${input_tiebrush_junc##*/}"
  base_name="${f%.*}"
fi

jscore_out="$outdir/${base_name}.jscore.txt"
tsstes_scores_out="$outdir/${base_name}.tsstes.scores.txt"
combined_out="$outdir/${base_name}.combined.ptf"

eval_junc_out="$outdir/${base_name}.eval.junctions.txt"
eval_tsstes_out="$outdir/${base_name}.eval.tsstes.txt"

# LightGBM threshold flags
declare -a score_flags=()
if [[ -n "$score_arg" ]]; then
  score_flags=(-s "$score_arg")
fi

# --- Small utilities for evaluation based on out/ files ---
filter_pos_from_jscore() {
  # input: jscore file (with predicted_label as last col)
  # output: lines with predicted_label==1, header dropped
  awk 'NR==1{next} ($NF==1 || $NF==1.0)' "$1"
}

filter_pos_from_tsstes_scores() {
  awk 'NR==1{next} ($NF==1 || $NF==1.0)' "$1"
}

build_annotation_splice_sites() {
  # outputs a 2-col file: chr <space/tab> pos
  local gtf="$1"
  local out="$2"
  local tmp_introns="$3"

  python3 "${helpers_dir}/gtf_to_intron_bed.py" "$gtf" "$tmp_introns"

  # intron bed is expected: chr start end ...
  awk 'BEGIN{OFS="\t"} {print $1,$2; print $1,$3}' "$tmp_introns" \
    | sort -k1,1 -k2,2n | uniq > "$out"
}

build_annotation_tsstes() {
  # outputs: chr \t pos \t TYPE  where TYPE in {TSS,CPAS}
  # Prefers 'transcript' features if present; else infers per transcript_id from exons.
  local gtf="$1"
  local out="$2"

  python3 - "$gtf" "$out" <<'PY'
import sys, re
gtf, outp = sys.argv[1], sys.argv[2]

attr_re = re.compile(r'(\S+)\s+"([^"]+)"')

def parse_attrs(s):
    d={}
    for m in attr_re.finditer(s):
        d[m.group(1)] = m.group(2)
    return d

# We will collect:
# - transcript features if present: tid -> (chr,strand,start,end)
# - else exon extrema per transcript_id: tid -> min_start, max_end, chr,strand
tx = {}
exon_ext = {}

with open(gtf, "r", encoding="utf-8", errors="replace") as f:
    for line in f:
        if not line or line[0] == "#":
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        chrom, source, feature, start, end, score, strand, frame, attrs = parts
        try:
            start_i = int(start)
            end_i = int(end)
        except:
            continue
        a = parse_attrs(attrs)
        tid = a.get("transcript_id") or a.get("transcriptId") or a.get("transcript")  # last one is rare
        if not tid:
            continue

        if feature == "transcript":
            tx[tid] = (chrom, strand, start_i, end_i)
        elif feature == "exon":
            v = exon_ext.get(tid)
            if v is None:
                exon_ext[tid] = [chrom, strand, start_i, end_i]
            else:
                # keep chrom/strand as first-seen; update bounds
                v[2] = min(v[2], start_i)
                v[3] = max(v[3], end_i)

# Choose transcript coords source
coords = tx if tx else exon_ext

# Emit unique positions
seen = set()
out = []
for tid, (chrom, strand, s, e) in coords.items():
    if strand == "+":
        tss = s
        cpas = e
    else:
        tss = e
        cpas = s
    k1 = (chrom, tss, "TSS")
    k2 = (chrom, cpas, "CPAS")
    if k1 not in seen:
        seen.add(k1); out.append(k1)
    if k2 not in seen:
        seen.add(k2); out.append(k2)

out.sort(key=lambda x: (x[0], x[1], x[2]))
with open(outp, "w") as w:
    for chrom, pos, typ in out:
        w.write(f"{chrom}\t{pos}\t{typ}\n")
PY
}

run_evaluation_from_outputs() {
  # Requires: out/<base>.jscore.txt and out/<base>.tsstes.scores.txt exist
  # Requires: -a provided
  local gtf="$1"
  local base="$2"

  local workdir
  workdir="$(mktemp -d "${outdir}/.eval.${base}.XXXX")"
  trap 'rm -rf "$workdir"' RETURN

  local ann_introns="$workdir/ann.introns.bed"
  local ann_ss="$workdir/ann.uniquess.bed"
  local ann_tsstes="$workdir/ann.tsstes.bed"

  local jpos="$workdir/${base}.jscore.pos.txt"
  local tspos="$workdir/${base}.tsstes.pos.txt"

  log "Eval: building annotation references from GTF..."
  build_annotation_splice_sites "$gtf" "$ann_ss" "$ann_introns"
  build_annotation_tsstes "$gtf" "$ann_tsstes"

  log "Eval: extracting predicted positives from out/ scores..."
  filter_pos_from_jscore "$jscore_out" > "$jpos" || true
  filter_pos_from_tsstes_scores "$tsstes_scores_out" > "$tspos" || true

  log "Eval: junctions -> ${eval_junc_out}"
  if [[ -s "$jpos" && -s "$ann_ss" ]]; then
    # evaluate_ptf_new.py prints to stdout; capture to file
    python3 "${helpers_dir}/evaluate_ptf_new.py" "$jpos" "$ann_ss" > "$eval_junc_out" || true
  else
    printf "No junction positives or empty annotation splice-site set.\n" > "$eval_junc_out"
  fi

  log "Eval: TSSTES -> ${eval_tsstes_out}"
  if [[ -s "$tspos" && -s "$ann_tsstes" ]]; then
    python3 "${helpers_dir}/evaluate_TSSTES_new.py" "$tspos" "$ann_tsstes" > "$eval_tsstes_out" || true
  else
    printf "No TSSTES positives or empty annotation TSSTES set.\n" > "$eval_tsstes_out"
  fi

  log "Eval complete."
}

# ---------------------------
# EVAL-ONLY MODE
# ---------------------------
if $eval_only; then
  log "Mode: eval-only"
  [[ -f "$jscore_out" ]] || die "Missing: $jscore_out (run full pipeline first, or use correct -b)"
  [[ -f "$tsstes_scores_out" ]] || die "Missing: $tsstes_scores_out (run full pipeline first, or use correct -b)"
  run_evaluation_from_outputs "$input_annotation" "$base_name"
  log "Done (eval-only)."
  exit 0
fi

# ---------------------------
# FULL PIPELINE MODE
# ---------------------------
log "Mode: full pipeline"
need_cmd bigWigToBedGraph

[[ -f "$input_tiebrush_junc" ]]   || die "Junction file not found: $input_tiebrush_junc"
[[ -f "$input_tiebrush_bigwig" ]] || die "BigWig file not found: $input_tiebrush_bigwig"

workdir="$(mktemp -d "${outdir}/.work.${base_name}.XXXX")"
cleanup() {
  if [[ -n "${KEEP_TEMP:-}" ]]; then
    log "KEEP_TEMP set; leaving workspace: $workdir"
  else
    rm -rf "$workdir"
  fi
}
trap cleanup EXIT

sorted_junc="$workdir/${base_name}.sorted.bed"
processed_junc="$workdir/${base_name}.jproc.txt"
processed_junc_bundle="$workdir/${base_name}.jbund.txt"
jpos_source="$workdir/${base_name}.jpos.txt"
tsstes_pos_tmp="$workdir/${base_name}.tsstes.pos.txt"
jpos_ptf_tmp="$workdir/${base_name}.jpos.ptf"
converted_bedgraph="$workdir/${base_name}.bw.bedGraph"
round1_processed_junc="$workdir/${base_name}.r1.jproc.txt"
round2_processed_bundles="$workdir/${base_name}.bund.txt"
round2_processed_bundles_w_metrics="$workdir/${base_name}.r2.metrics.txt"
round2_processed_bundles_w_metrics_ptf="$workdir/${base_name}.r2.metrics.ptf"
round2_processed_bundles_w_metrics_tsstes_ptf="$workdir/${base_name}.tsstes.ptf"

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

log "Step 4: LightGBM scoring (junctions) -> ${jscore_out}"
python3 "${helpers_dir}/LightGBM_no_normscale.py" \
  -i "$processed_junc_bundle" \
  -o "$jscore_out" \
  ${score_flags[@]+"${score_flags[@]}"}

log "Step 5: Filtering score-positive junctions (temp)..."
awk 'NR==1{next} ($NF==1 || $NF==1.0)' "$jscore_out" > "$jpos_source" || true

log "Step 7: Emit PTF (junctions) -> temp"
awk 'BEGIN{OFS="\t"} { print $1, $2, $5, $12 }' "$jpos_source" > "$jpos_ptf_tmp" || true

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

log "Step 13: Filter TSSTES score-positive -> temp"
awk 'NR==1{next} ($NF==1 || $NF==1.0)' "$tsstes_scores_out" > "$tsstes_pos_tmp" || true

log "Step 15: Combine ptfs -> ${combined_out}"
"${helpers_dir}/combine_ptfs.sh" \
  "$jpos_source" \
  "$tsstes_pos_tmp" \
  > "$combined_out"

log "Final outputs:"
ls -lh "$jscore_out" "$tsstes_scores_out" "$combined_out" || true

# --- NEW: evaluation that uses generated outputs in out/ ---
if $anno_present; then
  log "Evaluation requested (-a): evaluating based on generated out/ files..."
  run_evaluation_from_outputs "$input_annotation" "$base_name"
fi

log "spliceCOV complete!"
