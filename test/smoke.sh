#!/usr/bin/env bash
set -euo pipefail

# This script fabricates tiny junctions/coverage/gtf inputs, then runs a minimal pipeline.
# Requirements (provided by CI steps):
#   - bedGraphToBigWig  (downloaded per OS)
#   - bigWigToBedGraph  (already installed by CI steps)
#   - bash, coreutils, awk

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="${TMPDIR:-/tmp}/splicecov-mini-$RANDOM"
mkdir -p "$TMP"

echo "[smoke] workdir: $TMP"

# 1) Minimal chrom sizes (1 Mb chr1)
cat > "$TMP/chrom.sizes" <<EOF
chr1	1000000
EOF

# 2) Minimal bedGraph coverage (two simple blocks)
cat > "$TMP/coverage.bedGraph" <<EOF
chr1	1000	2000	5
chr1	3000	3500	12
EOF

# 3) Build tiny bigWig from bedGraph
: "${BEDGRAPH_TO_BIGWIG:=bedGraphToBigWig}"
$BEDGRAPH_TO_BIGWIG "$TMP/coverage.bedGraph" "$TMP/chrom.sizes" "$TMP/coverage.bigWig"

# 4) Tiny junctions file (6-column BED-ish)
cat > "$TMP/junctions.bed" <<'EOF'
chr1	1200	1500	JUNC00000001	3	+
chr1	3200	3400	JUNC00000002	7	-
EOF

# 5) Tiny GTF (optional; two simple exons)
cat > "$TMP/annotation.gtf" <<'EOF'
chr1	test	source	exon	1100	1300	.	+	.	gene_id "g1"; transcript_id "t1";
chr1	test	source	exon	1400	1600	.	+	.	gene_id "g1"; transcript_id "t1";
chr1	test	source	exon	3100	3300	.	-	.	gene_id "g2"; transcript_id "t2";
chr1	test	source	exon	3350	3450	.	-	.	gene_id "g2"; transcript_id "t2";
chr1	test	source	transcript  3350	3450	.	-	.	gene_id "g2"; transcript_id "t2";
chr1	test	source	transcript	3350	3450	.	-	.	gene_id "g2"; transcript_id "t2";
EOF

echo "[smoke] Inputs created:"
ls -lh "$TMP"

# 6) Where the installed launcher lives (in CI we set PREFIX to a temp dir)
BIN_DIR="${BIN_DIR:-$HOME/.local/bin}"
LAUNCHER="${BIN_DIR}/splicecov"

echo "[smoke] Using launcher: $LAUNCHER"
"$LAUNCHER" -h >/dev/null || true

# Minimal run (junction + coverage)
echo "[smoke] run: minimal (junction + coverage)"
"$LAUNCHER" -j "$TMP/junctions.bed" -c "$TMP/coverage.bigWig"

# With annotation to exercise eval path if present
echo "[smoke] run: with annotation"
"$LAUNCHER" -j "$TMP/junctions.bed" -c "$TMP/coverage.bigWig" -a "$TMP/annotation.gtf"

echo "[smoke] OK"
