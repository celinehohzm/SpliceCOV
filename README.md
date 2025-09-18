# SpliceCOV

SpliceCOV uses a **LightGBM-based machine learning model** to score and predict splice sites and TSS/TESs from RNA-seq coverage evidence.

By providing these predictions, SpliceCOV improves **precision while maintaining sensitivity**, especially when running StringTie with the --ptf option.

It can also evaluate predictions against a genome annotation (GTF) if provided.

---

## Getting Started

```bash
git clone https://github.com/celinehohzm/SpliceCOV.git
cd SpliceCOV/scripts

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

### Requirements

- `bash` (version 4 or newer)
- Standard Linux tools: `awk`, `sort`, `comm` (comm only if using `-a`)
- UCSC `bigWigToBedGraph` (from UCSC Genome Browser utilities)
  - Linux/macOS binaries: `http://hgdownload.soe.ucsc.edu/admin/exe/`
  - Put the binary in your `PATH` (e.g., `/usr/local/bin/`).
  - macOS (Homebrew): `brew install ucsc-genome-browser`
- Python (version 3.8 or newer) with:
  - `lightgbm`
  - `numpy`
  - `pandas`
- **Optional (recommended):** `GNU time` for detailed runtime/memory stats  
  - macOS: `brew install gnu-time` (available as `gtime`)

### Python setup with conda

```bash
conda create -n splicecov python=3.10 # Python version 3.10
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

**Full CLI:**
```
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
```

---

**How basename -b is chosen**

SpliceCOV writes all outputs to out/ and prefixes them with a basename:

1. If you pass -b <name>, SpliceCOV uses exactly <name> as the basename.
   
- Allowed: letters, numbers, underscores, dashes (no slashes/spaces).
- Example: -b `gtex_v8_spleen`

2. If no `-b` is given, SpliceCOV derives the basename from the junctions file given to `-j`:

- Example: `-j /data/juncs/sample.txt`. then, `$basename` = `sample`

All output files will be named like `out/${basename}.*`.


---
**What does -s (score threshold) do?**

- `-s` sets the probability threshold used when converting LightGBM scores to positive predictions.

- It is forwarded to:

  - Step 4: junction scoring (LightGBM_no_normscale.py)

  - Step 12: TSS/CPAS scoring (LightGBM_tss.py)

- Valid range is 0–1. If you omit -s, the Python scripts use their default (0.4).


**Examples:**
```
# Stricter calling (fewer positives)
./splicecov.sh -j sample.tiebrush_junctions.txt -c sample.coverage.bigWig -s 0.8

# More permissive calling (more positives)
./splicecov.sh -j sample.tiebrush_junctions.txt -c sample.coverage.bigWig -s 0.25
```

---
## Inputs (recommended to generate from TieBrush & TieCov)

SpliceCOV expects:
- A TieBrush junctions file (-j) — aggregated splice junctions
- A BigWig coverage file (-c) — aggregated per-base read coverage
These typically come from the TieBrush and TieCov tools (https://github.com/alevar/tiebrush).



Junctions file format (columns in order: chr, junction_start, junction_end, junction_name, number_of_samples, strand):
```
track name=junctions
chr1    10620   21193   JUNC00000001    1       +
chr1    11348   11410   JUNC00000002    1       -
chr1    11410   17721   JUNC00000003    1       -
chr1    11439   21424   JUNC00000004    1       -
chr1    11671   12009   JUNC00000005    12      +
```


Bedgraph file format generated from bigwig (colummns in order: chr, position_start, position_end, coverage):
```
chr1    10535   10538   1
chr1    10538   10540   2
chr1    10540   10542   4
chr1    10542   10545   7
chr1    10545   10560   8
```

---
## Output 
All outputs are written to the `out/` folder.

**Core junction outputs**

- `${basename}`.sorted.bed — sorted junctions
- `${basename}`.jproc.txt — processed junctions
- `${basename}`.jbund.txt — junctions + coverage
- `${basename}`.jscore.txt — LightGBM scores (junctions)
- `${basename}`.jpos.txt — score-positive junctions
- `${basename}`.jpos.ptf — PTF for score-positive junctions

**Round-2 / TSSTES outputs**

- `${basename}`.bw.bedGraph — BigWig converted to BedGraph
- `${basename}`.r2.bund.txt — round-2 bundles
- `${basename}`.r2.metrics.txt — TSSTES metrics
- `${basename}`.r2.metrics.ptf — PTF after metrics
- `${basename}`.r2.tsstes.ptf — TSS/CPAS only
- `${basename}`.r2.tsstes.scores.txt — LightGBM scores (TSS/CPAS)
- `${basename}`.r2.tsstes.pos.txt — score-positive TSS/CPAS
- `${basename}`.r2.tsstes.pos.eval.txt — evaluation vs annotation (only with -a)

**Final combined output**

- `${basename}`.combined.ptf — combined score-positive splice-sites + TSS/TESs


---
## Running Stringtie with SpliceCOV's output

After generating the `${basename}.combined.ptf` file, you can run Stringtie with the `-ptf` option. 

Ideally it will boost precision while retaining sensitivity when running with SpliceCOV's `-ptf` predictions. 

Check out Stringtie's manual here: https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual

`stringtie sample.bam -ptf out/<basename>_combined.ptf -o sample.stringtie.gtf `

--- 
## Tips & Troubleshooting

- If you see “`bigWigToBedGraph: command not found`”, install the UCSC utilities and ensure they’re in your PATH.

- On macOS, GNU time appears as gtime; SpliceCOV auto-detects it for better runtime/memory summaries.

- Use `-b` to keep runs tidy and identifiable (e.g., -b `gtex_v8_brain_cortex`).

- Tune `-s` to adjust precision/recall in your calls; higher thresholds are stricter.

---
**Questions or issues?** Please open a GitHub issue with your command line, environment (OS, Python version, tool versions), and a short log excerpt.



