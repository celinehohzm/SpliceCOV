# SpliceCOV

SpliceCOV is a command-line pipeline that processes RNA-seq coverage and junction data to identify and score **splice junctions** and **transcription start/termination events (TSS/CPAS)**.  
It can also evaluate predictions against a genome annotation (GTF) if provided.

---

## Getting Started

```bash
git clone https://github.com/celinehohzm/SpliceCOV.git
cd SpliceCOV

# Minimal run (junction + coverage only)
./splicecov.sh -j sample.tiebrush_junctions.txt -c sample.coverage.bigWig

# Run with annotation (adds evaluation steps)
./splicecov.sh -j sample.tiebrush_junctions.txt -c sample.coverage.bigWig -a gencode.v43.annotation.gtf

# Resume from step 8
./splicecov.sh -j sample.tiebrush_junctions.txt -c sample.coverage.bigWig -n 8

# Show help
./splicecov.sh -h
```

---
## Installation

SpliceCOV is written in Bash with helper scripts in Perl and Python.

**Requirements**

- `bash` (≥4)
- Standard Linux tools: `awk`, `sort`, `comm` (comm only if using -a)
- `bigWigToBedGraph` (from UCSC Genome Browser utilities)
- Python ≥3.8 with:
  - `lightgbm`
  - `numpy`
  - `pandas`

Install Python deps with conda:
```conda create -n splicecov python=3.10
conda activate splicecov
pip install lightgbm numpy pandas
```

---
## Usage
**Junction-only mode**

If you have TieBrush junctions and a BigWig coverage file, run:

`./splicecov.sh -j sample.tiebrush_junctions.txt -c sample.coverage.bigWig`

**Annotation mode**

If you also provide a **GTF annotation**, SpliceCOV will build reference introns and unique splice sites, and evaluate predictions:

`./splicecov.sh -j sample.tiebrush_junctions.txt -c sample.coverage.bigWig -a gencode.v43.annotation.gtf`

**Resuming with `-n`**

SpliceCOV runs as numbered steps (1 → 15).
If interrupted, resume from a later step:

`./splicecov.sh -j sample.tiebrush_junctions.txt -c sample.coverage.bigWig -n 8`

---
## Output

All outputs are written to the `out/` folder.
Assume your input junction file is `sample.txt` → `<b> = sample`.

Core junction outputs

<b>.sorted.bed — sorted junctions

<b>.jproc.txt — processed junctions

<b>.jbund.txt — junctions + coverage

<b>.jscore.txt — LightGBM scores (junctions)

<b>.jpos.txt — score-positive junctions

<b>.jpos.ptf — PTF for score-positive junctions

Round-2 / TSSTES outputs

<b>.bw.bedGraph — BigWig converted to BedGraph

<b>.r2.bund.txt — round-2 bundles

<b>.r2.metrics.txt — TSSTES metrics

<b>.r2.metrics.ptf — PTF after metrics

<b>.r2.tsstes.ptf — TSS/CPAS only

<b>.r2.tsstes.scores.txt — LightGBM scores (TSS/CPAS)

<b>.r2.tsstes.pos.txt — score-positive TSS/CPAS

<b>.r2.tsstes.pos.eval.txt — evaluation vs annotation (only with -a)

Final combined output

<b>.combo.ptf.tsv — combined score-positive junctions + TSS/CPAS
