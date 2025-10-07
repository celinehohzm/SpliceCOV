# SpliceCOV: Accurate Identification of Transcriptional Features from RNA Sequencing

SpliceCOV is a segmentation pipeline for **detecting splice sites and transcription start/end sites (TSS/TES) from RNA-Seq coverage**. SpliceCOV follows a two-step approach: it first employs a fast heuristic to identify candidate changepoints, and then applies a LightGBM decision-tree classifier to filter and refine these candidates, generating high-confidence transcriptional boundaries. 

The required inputs for SpliceCOV are the junctions (in `bed` format) file and coverage (in `bigwig` format) file. It is recommended to generate these using Tiebrush/TieCOV (https://github.com/alevar/tiebrush).

When SpliceCOV's predictions are used to guide transcriptome assembly with StringTie's  `--ptf` option (https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual),  the resulting transcripts exhibit increased precision without sacrificing sensitivity. As the number of samples per tissue increases, SpliceCOV continues to improve its accuracy in detecting transcripts, highlighting its scalability. 

---

## Installation

You can build from source by:
```bash
git clone https://github.com/celinehohzm/SpliceCOV.git
cd SpliceCOV
conda env create -f /path/to/SpliceCOV/environment.yml
conda activate splicecov
make PREFIX="$HOME/.local" release

echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc   # Linux/bash
# or
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.zshrc    # macOS/zsh
```

Sample commands: 

```bash
# Minimal run (junction + coverage only)
splicecov.sh -j sample.tiebrush_junctions.txt -c sample.coverage.bigWig

# Run with annotation (adds evaluation steps)
splicecov.sh -j sample.tiebrush_junctions.txt -c sample.coverage.bigWig -a gencode.v43.annotation.gtf

# Show help
splicecov.sh -h
```

---
## Installation Requirements

SpliceCOV is written in Bash with helper scripts in Perl and Python. The softwares needed are:

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

---
## Inputs (recommended to generate with TieBrush & TieCov)

SpliceCOV expects the following:

- TieBrush junctions file (`-j`) — aggregated splice junctions
- BigWig coverage file (`-c`) — aggregated per-base read coverage

Both of these files are typically produced by the TieBrush and TieCov tools (https://github.com/alevar/tiebrush).



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
## Usage
**Junction-only mode**

If you have TieBrush junctions and a BigWig coverage file, run:

`./splicecov.sh -j sample.tiebrush_junctions.txt -c sample.coverage.bigWig`

**Annotation mode**

If you also provide a **GTF annotation**, SpliceCOV will build reference introns and unique splice sites, and evaluate predictions:

`./splicecov.sh -j sample.tiebrush_junctions.txt -c sample.coverage.bigWig -a gencode.v43.annotation.gtf`

**Full CLI:**
```
Usage: splicecov.sh -j <input_tiebrush_junc> -c <input_tiebrush_bigwig> [-a <annotation_gtf>] [-b <basename>] [-s <threshold>]

Required:
  -j <file> : input TieBrush junction file 
  -c <file> : input coverage BigWig file 

Optional:
  -a <file> : input annotation (GTF). If provided, annotation-dependent steps
              (building introns/unique splice sites and evaluation) will run.
  -b <str>  : basename to use for ALL output files; overrides the default from -j.
  -s <num>  : LightGBM scoring threshold [0,1], default is 0.4.

Outputs (only these remain in out/):
  <basename>.jscore.txt
  <basename>.tsstes.scores.txt
  <basename>.combined.ptf
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
## Output 
All outputs are written to the `out/` folder.

**Core outputs**

- `${basename}`.jscore.txt — LightGBM scores (junctions)
- `${basename}`.tsstes.scores.txt — LightGBM scores (TSS/CPAS)
- `${basename}`.combined.ptf — combined score-positive splice-sites + TSS/TESs


---
## Running Stringtie with SpliceCOV's output

After generating the `${basename}.combined.ptf` file, you can run Stringtie with the `-ptf` option. 

Check out Stringtie's manual here: https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual

Sample command:
`stringtie sample.bam --ptf out/<basename>_combined.ptf -o sample.stringtie.gtf `

As shown in the figure below, SpliceCOV-guided assemblies consistently achieve higher precision without compromising sensitivity, as measured by the number of reference‐matching transcripts. Moreover, the precision gain increases with tissue cohort size, likely reflecting the cleaner coverage signals in larger sample sets.

<img width="637" height="794" alt="image" src="https://github.com/user-attachments/assets/d186eed9-be27-40d9-b6d4-b610ce246eaf" />

 

--- 
## Tips & Troubleshooting

- If you see “`bigWigToBedGraph: command not found`”, install the UCSC utilities and ensure they’re in your PATH.

- On macOS, GNU time appears as gtime; SpliceCOV auto-detects it for better runtime/memory summaries.

- Use `-b` to keep runs tidy and identifiable (e.g., -b `gtex_v8_brain_cortex`).

- Tune `-s` to adjust precision/recall in your calls; higher thresholds are stricter.

---
**Questions or issues?** Please open a GitHub issue with your command line, environment (OS, Python version, tool versions), and a short log excerpt.



