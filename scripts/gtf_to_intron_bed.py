import sys

def parse_gtf(gtf_file):
    """Parse GTF file and extract intron junctions."""
    transcripts = {}
    
    with open(gtf_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if parts[2] == "exon":
                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                transcript_id = None
                
                # Extract transcript_id
                attributes = parts[8].split(';')
                for attr in attributes:
                    attr = attr.strip()
                    if attr.startswith('transcript_id'):
                        transcript_id = attr.split('"')[1]
                        break
                
                if transcript_id:
                    if transcript_id not in transcripts:
                        transcripts[transcript_id] = []
                    transcripts[transcript_id].append((chrom, start, end, strand))
    
    return transcripts

def extract_introns(transcripts):
    """Extract introns from transcript exon data."""
    introns = []
    
    for transcript_id, exons in transcripts.items():
        # Sort exons by start position
        exons.sort(key=lambda x: x[1])
        
        for i in range(len(exons) - 1):
            chrom, start1, end1, strand = exons[i]
            _, start2, _, _ = exons[i + 1]
            
            intron_start = end1 
            intron_end = start2 
            introns.append((chrom, intron_start, intron_end, strand, transcript_id))
    
    return introns

def write_bed(introns, bed_file):
    """Write introns to BED file."""
    with open(bed_file, 'w') as file:
        for intron in introns:
            chrom, start, end, strand, transcript_id = intron
            file.write(f"{chrom}\t{start}\t{end}\t{transcript_id}\t0\t{strand}\n")

def main(gtf_file, bed_file):
    transcripts = parse_gtf(gtf_file)
    introns = extract_introns(transcripts)
    write_bed(introns, bed_file)
    print(f"Extracted {len(introns)} introns and saved to {bed_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_introns.py <input.gtf> <output.bed>")
        sys.exit(1)
    
    gtf_file = sys.argv[1]
    bed_file = sys.argv[2]
    
    main(gtf_file, bed_file)
