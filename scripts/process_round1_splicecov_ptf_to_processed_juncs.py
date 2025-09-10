import sys
from collections import defaultdict

def process_ptf_file(input_file, additional_bed_file, junctions_bed_file):
    """
    Processes a .ptf file to generate a junctions.bed file.
    Matches JSTART and JEND rows based on their junction name.
    If a lonely row is found, tries to find a match in the additional BED file and includes it in the output.
    """
    # Dictionaries to store JSTART and JEND rows by junction name
    jstart_rows = {}
    jend_rows = {}
    additional_bed_data = {}

    # Read the additional BED file into a dictionary
    with open(additional_bed_file, 'r') as bed_file:
        for line in bed_file:
            # Skip comments or empty lines
            if not line.strip() or line.startswith("#"):
                continue

            # Parse the line
            columns = line.strip().split("\t")
            if len(columns) < 4:
                # Ensure there are at least 4 columns to avoid IndexError
                continue
            junction_name = columns[3]  # Column 4: Junction name
            additional_bed_data[junction_name] = line.strip()

    # Read the input .ptf file
    with open(input_file, 'r') as infile:
        for line_num, line in enumerate(infile, 1):
            # Skip comments or empty lines
            if not line.strip() or line.startswith("#"):
                continue

            # Parse the line
            columns = line.strip().split("\t")
            if len(columns) < 10:
                # Skip lines that do not have at least 10 columns
                print(f"Warning: Line {line_num} in {input_file} has fewer than 10 columns. Skipping.")
                continue

            # Truncate to first 10 columns to avoid unpacking errors
            columns = columns[:10]
            try:
                chrom, pos, junction_name, score, strand, percsame, _, _, _, flag = columns
                pos = int(pos)
            except ValueError as ve:
                print(f"Error parsing line {line_num} in {input_file}: {ve}. Skipping.")
                continue

            # Separate JSTART and JEND
            if flag == "JSTART":
                jstart_rows[junction_name] = columns
            elif flag == "JEND":
                jend_rows[junction_name] = columns

    # Open output file
    with open(junctions_bed_file, 'w') as junctions_bed:
        # Write the header for the junctions.bed file
        junctions_bed.write(
            '# track name=junctions type=bedDetail description="percsame-percdifferent-percleft-percright"\n'
        )

        # Match JSTART and JEND by junction name
        matched_junctions = set()
        for junction_name in jstart_rows:
            if junction_name in jend_rows:
                # Extract relevant fields for matched rows
                jstart = jstart_rows[junction_name]
                jend = jend_rows[junction_name]
                chrom = jstart[0]
                start = jstart[1]
                end = jend[1]
                score = jstart[3]
                strand = jstart[4]
                percsame = jstart[5]

                # Write the matched row to the junctions.bed file
                junctions_bed.write(f"{chrom}\t{start}\t{end}\t{junction_name}\t{score}\t{strand}\t{percsame}\n")
                matched_junctions.add(junction_name)

        # Include rows from the additional BED file for unmatched junctions
        for junction_name, row in jstart_rows.items():
            if junction_name not in matched_junctions and junction_name in additional_bed_data:
                junctions_bed.write(additional_bed_data[junction_name] + "\n")

        for junction_name, row in jend_rows.items():
            if junction_name not in matched_junctions and junction_name in additional_bed_data:
                junctions_bed.write(additional_bed_data[junction_name] + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python process_ptf.py <input_ptf_file> <additional_bed_file> <junctions_bed_file>")
        sys.exit(1)

    input_ptf_file = sys.argv[1]
    additional_bed_file = sys.argv[2]
    junctions_bed_file = sys.argv[3]

    process_ptf_file(input_ptf_file, additional_bed_file, junctions_bed_file)
