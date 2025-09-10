SpliceCOV

SpliceCOV is a command-line pipeline that scores splice junctions and transcription start/termination events (TSS/CPAS) from RNA-seq coverage and (optionally) evaluates predictions against a reference annotation.

If you’re new to this:

A splice junction is where exons join after splicing.

A BigWig is a compact file showing read coverage along the genome.

SpliceCOV combines junction coordinates with coverage to identify well-supported junctions and TSS/CPAS.

Quick start
# Minimal run (no annotation)
./splicecov.sh -j path/to/tiebrush_junctions.txt -c path/to/coverage.bigWig

# With annotation (adds evaluation steps)
./splicecov.sh -j path/to/tiebrush_junctions.txt -c path/to/coverage.bigWig -a path/to/annotation.gtf


All outputs land in out/. The final combined PTF is:

out/<basename>.combo.ptf.tsv


where <basename> is your junction filename without the .txt suffix.

Requirements

Bash, awk, sort, comm (comm used only when -a is provided)

bigWigToBedGraph (UCSC tools) for BigWig → BedGraph conversion

Python 3 with packages used by helper scripts (e.g., LightGBM)

(Optional) /usr/bin/time or gtime for a run summary
