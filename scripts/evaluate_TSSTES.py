#!/usr/bin/env python3

import sys
import os
import bisect

# Reference file path (fixed)
ref_file_path = "/home/choh1/spliceCOV/data/chess.primary.unique_tsstes.bed"

def read_reference_file(ref_file_path):
    """
    Reads the reference BED file and stores positions in sorted lists per chromosome.
    """
    ref_dict = {}
    try:
        with open(ref_file_path, 'r') as ref_file:
            for line in ref_file:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                cols = line.split()
                if len(cols) >= 2:
                    chrom = cols[0]
                    try:
                        pos = int(cols[1])
                    except ValueError:
                        print(f"Warning: Non-integer position '{cols[1]}' in reference file. Skipping line.")
                        continue
                    if chrom not in ref_dict:
                        ref_dict[chrom] = []
                    ref_dict[chrom].append(pos)
        # Sort the positions for each chromosome
        for chrom in ref_dict:
            ref_dict[chrom] = sorted(ref_dict[chrom])
        return ref_dict
    except FileNotFoundError:
        print(f"Error: Reference file '{ref_file_path}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while reading the reference file: {e}")
        sys.exit(1)

def find_matches(target_rows, ref_dict, window=50):
    """
    Finds matching and non-matching rows from the target file.
    """
    matching_rows = []
    non_matching_rows = []
    matched_ref_positions = set()

    for row in target_rows:
        chrom = row['chrom']
        pos = row['pos']
        if chrom not in ref_dict:
            non_matching_rows.append(row)
            continue
        ref_positions = ref_dict[chrom]
        left = pos - window
        right = pos + window
        left_idx = bisect.bisect_left(ref_positions, left)
        right_idx = bisect.bisect_right(ref_positions, right)
        if left_idx < right_idx:
            matching_rows.append(row)
            for ref_pos in ref_positions[left_idx:right_idx]:
                matched_ref_positions.add((chrom, ref_pos))
        else:
            non_matching_rows.append(row)

    return matching_rows, non_matching_rows, matched_ref_positions

def add_label_column(input_file_path, output_file_path, ref_dict, window=50):
    """
    Creates a new file with an additional column: "1" for matching rows, "0" for non-matching rows.
    """
    try:
        with open(input_file_path, 'r') as infile, open(output_file_path, 'w') as outfile:
            for line in infile:
                line = line.strip()
                if not line or line.startswith('#'):
                    outfile.write(line + '\n')
                    continue
                cols = line.split()
                if len(cols) < 2:
                    continue
                chrom = cols[0]
                try:
                    pos = int(cols[1])
                except ValueError:
                    outfile.write(line + '\t0\n')
                    continue
                if chrom not in ref_dict:
                    outfile.write(line + '\t0\n')
                    continue
                ref_positions = ref_dict[chrom]
                left = pos - window
                right = pos + window
                left_idx = bisect.bisect_left(ref_positions, left)
                right_idx = bisect.bisect_right(ref_positions, right)
                if left_idx < right_idx:
                    outfile.write(line + '\t1\n')
                else:
                    outfile.write(line + '\t0\n')
        print(f"File with labeled rows written to: {output_file_path}")
    except Exception as e:
        print(f"An error occurred while adding labels to rows: {e}")
        sys.exit(1)

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <target_file>")
        sys.exit(1)

    target_file_path = sys.argv[1]
    target_file_basename = os.path.basename(target_file_path)

    # Original output files
    matching_output_file = f"{target_file_basename}_matching_rows.txt"
    non_matching_output_file = f"{target_file_basename}_non_matching_clusters.txt"
    ref_only_output_file = f"{target_file_basename}_ref_only_rows.txt"

    # New output file
    labeled_output_file = f"{target_file_basename}_with_labels.txt"

    # Read reference file
    print("Reading reference BED file...")
    ref_dict = read_reference_file(ref_file_path)

    # Add labeled column to the input file
    print(f"Processing target file '{target_file_path}' to add labeled column...")
    add_label_column(target_file_path, labeled_output_file, ref_dict)

    # Process original outputs
    print("Generating original outputs...")
    target_rows = []
    try:
        with open(target_file_path, 'r') as target_file:
            for line in target_file:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                cols = line.split()
                if len(cols) < 2:
                    continue
                chrom = cols[0]
                try:
                    pos = int(cols[1])
                except ValueError:
                    continue
                target_rows.append({'chrom': chrom, 'pos': pos, 'line': line})
    except Exception as e:
        print(f"An error occurred while reading the target file: {e}")
        sys.exit(1)

    matching_rows, non_matching_rows, matched_ref_positions = find_matches(target_rows, ref_dict, window=50)
    print(f"Found {len(matching_rows)} matching rows and {len(non_matching_rows)} non-matching rows.")

    try:
        with open(matching_output_file, 'w') as match_file:
            for row in matching_rows:
                match_file.write(row['line'] + '\n')
        print(f"Matching rows written to: {matching_output_file}")
    except Exception as e:
        print(f"An error occurred while saving matching rows: {e}")
        sys.exit(1)

    try:
        with open(ref_only_output_file, 'w') as ref_only_file:
            unmatched_count = 0
            for chrom in ref_dict:
                ref_positions = ref_dict[chrom]
                for pos in ref_positions:
                    if (chrom, pos) not in matched_ref_positions:
                        ref_only_file.write(f"{chrom}\t{pos}\n")
                        unmatched_count += 1
        print(f"Reference-only rows written to: {ref_only_output_file}")
        print(f"Number of reference rows with matches: {len(matched_ref_positions)}")
        print(f"Number of reference rows without matches: {unmatched_count}")

        # Calculate Precision and Sensitivity
        TP = len(matching_rows)
        FP = len(non_matching_rows)
        FN = unmatched_count

        precision = TP / (TP + FP) if (TP + FP) > 0 else 0
        sensitivity = TP / (TP + FN) if (TP + FN) > 0 else 0

        print(f"Precision: {precision:.4f}")
        print(f"Sensitivity: {sensitivity:.4f}")
    except Exception as e:
        print(f"An error occurred while saving reference-only rows: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
