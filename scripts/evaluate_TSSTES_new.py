#!/usr/bin/env python3

import sys
import os
import bisect

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
                    ref_dict.setdefault(chrom, []).append(pos)
        # Sort the positions for each chromosome
        for chrom in ref_dict:
            ref_dict[chrom].sort()
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
        chrom, pos = row['chrom'], row['pos']
        if chrom not in ref_dict:
            non_matching_rows.append(row)
            continue
        ref_positions = ref_dict[chrom]
        left, right = pos - window, pos + window
        l_idx = bisect.bisect_left(ref_positions, left)
        r_idx = bisect.bisect_right(ref_positions, right)
        if l_idx < r_idx:
            matching_rows.append(row)
            for ref_pos in ref_positions[l_idx:r_idx]:
                matched_ref_positions.add((chrom, ref_pos))
        else:
            non_matching_rows.append(row)

    return matching_rows, non_matching_rows, matched_ref_positions

def add_label_column(input_file_path, ref_dict, output_file_path, window=50):
    """
    Writes out the input file with an extra final column: "1" if within `window` of a ref position, else "0".
    """
    try:
        with open(input_file_path, 'r') as infile, open(output_file_path, 'w') as outfile:
            for line in infile:
                line = line.rstrip('\n')
                if not line or line.startswith('#'):
                    outfile.write(line + '\n')
                    continue
                cols = line.split()
                if len(cols) < 2:
                    outfile.write(line + '\t0\n')
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
                left, right = pos - window, pos + window
                l_idx = bisect.bisect_left(ref_positions, left)
                r_idx = bisect.bisect_right(ref_positions, right)
                label = "1" if l_idx < r_idx else "0"
                outfile.write(f"{line}\t{label}\n")
        print(f"Labeled file written to: {output_file_path}")
    except Exception as e:
        print(f"An error occurred while adding labels: {e}")
        sys.exit(1)

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {os.path.basename(sys.argv[0])} <target_file> <reference_bed>")
        sys.exit(1)

    target_file_path = sys.argv[1]
    ref_file_path    = sys.argv[2]
    basename         = os.path.basename(target_file_path)

    # Output filenames
    matching_output_file   = f"{basename}_matching_rows.txt"
    non_matching_output    = f"{basename}_non_matching_clusters.txt"
    ref_only_output_file   = f"{basename}_ref_only_rows.txt"
    labeled_output_file    = f"{basename}_with_labels.txt"

    # Read reference
    print("Reading reference BED file...")
    ref_dict = read_reference_file(ref_file_path)

    # Add label column
    print("Adding label column to target file...")
    add_label_column(target_file_path, ref_dict, labeled_output_file)

    # Load target rows for matching logic
    target_rows = []
    try:
        with open(target_file_path, 'r') as tf:
            for line in tf:
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
        print(f"Error reading target file: {e}")
        sys.exit(1)

    matching_rows, non_matching_rows, matched_ref_positions = find_matches(target_rows, ref_dict)

    print(f"Found {len(matching_rows)} matching and {len(non_matching_rows)} non-matching rows.")

    # Write matches
    try:
        with open(matching_output_file, 'w') as mf:
            for r in matching_rows:
                mf.write(r['line'] + '\n')
        print(f"Matching rows → {matching_output_file}")
    except Exception as e:
        print(f"Error writing matching rows: {e}")
        sys.exit(1)

    # Write reference-only (false negatives)
    try:
        fn_count = 0
        with open(ref_only_output_file, 'w') as rf:
            for chrom, positions in ref_dict.items():
                for pos in positions:
                    if (chrom, pos) not in matched_ref_positions:
                        rf.write(f"{chrom}\t{pos}\n")
                        fn_count += 1
        print(f"Reference-only rows → {ref_only_output_file}")
    except Exception as e:
        print(f"Error writing reference-only rows: {e}")
        sys.exit(1)

    # Compute and print metrics
    TP = len(matching_rows)
    FP = len(non_matching_rows)
    FN = fn_count

    precision   = TP / (TP + FP) if (TP + FP) > 0 else 0
    sensitivity = TP / (TP + FN) if (TP + FN) > 0 else 0

    print(f"Precision:   {precision:.4f}")
    print(f"Sensitivity: {sensitivity:.4f}")

if __name__ == "__main__":
    main()
