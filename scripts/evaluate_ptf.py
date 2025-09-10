#!/usr/bin/env python3

import sys
import os

# Reference file path (fixed)
# ref_file_path = "/home/choh1/spliceCOV/data/liver/REF_liver_tiebrush_in_chess.ss.bed"
# ref_file_path = "/home/choh1/spliceCOV/data/liver/REF_liver_tiebrush_in_chess.ss.test.bed"

ref_file_path = "/home/choh1/spliceCOV/data/bladder/REF_bladder_tiebrush_in_chess.ss.bed"
# ref_file_path = "/home/choh1/spliceCOV/data/bladder/REF_bladder_tiebrush_in_chess.ss.test.bed"

# ref_file_path = "/home/choh1/spliceCOV/data/small_intestine/REF_small_intestine_tiebrush_in_chess.ss.bed"
# ref_file_path = "/home/choh1/spliceCOV/data/small_intestine/REF_small_intestine_tiebrush_in_chess.ss.test.bed"
# ref_file_path = "/home/choh1/spliceCOV/data/small_intestine/REF_small_intestine_tiebrush_in_refseq.ss.bed"

# ref_file_path = "/home/choh1/spliceCOV/data/bone_marrow/REF_bone_marrow_tiebrush_in_chess.ss.bed"
# ref_file_path = "/home/choh1/spliceCOV/data/spleen/REF_spleen_tiebrush_in_chess.ss.bed"
# ref_file_path = "/home/choh1/spliceCOV/data/blood/REF_blood_tiebrush_in_chess.ss.bed"
# ref_file_path = '/home/choh1/spliceCOV/data/kidney/REF_kidney_tiebrush_in_chess.ss.bed'
# ref_file_path = "/home/choh1/spliceCOV/data/small_intestine/REF_small_intestine_tiebrush_in_gencode.ss.bed"

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <target_file>")
        sys.exit(1)

    target_file_path = sys.argv[1]
    target_file_basename = os.path.basename(target_file_path)

    # Define output file names
    matching_output_file = f"{target_file_basename}_matching_rows.txt"
    non_matching_output_file = f"{target_file_basename}_non_matching_rows.txt"
    ref_only_output_file = f"{target_file_basename}_ref_only_rows.txt"
    labeled_output_file = f"{target_file_basename}_with_labels.txt"

    # Read reference file and build a dict of positions to lines
    ref_positions = {}
    try:
        with open(ref_file_path, 'r') as ref_file:
            for line in ref_file:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                cols = line.split()
                if len(cols) >= 2:
                    chrom = cols[0]
                    pos = cols[1]
                    ref_positions[(chrom, pos)] = line
    except FileNotFoundError:
        print(f"Error: Reference file '{ref_file_path}' not found.")
        sys.exit(1)

    # Initialize counters and containers
    total_rows = 0
    matching_rows = {}
    non_matching_rows = {}
    target_positions = {}
    confidence_score_present = False

    try:
        with open(target_file_path, 'r') as target_file, open(labeled_output_file, 'w') as labeled_file:
            for line in target_file:
                line = line.strip()
                if not line or line.startswith('#'):
                    labeled_file.write(line + '\n')
                    continue
                cols = line.split()
                if len(cols) >= 2:
                    total_rows += 1
                    chrom = cols[0]
                    pos = cols[1]
                    key = (chrom, pos)

                    # Check dynamically for the confidence score column
                    if not confidence_score_present and len(cols) > 6:
                        try:
                            float(cols[6])
                            confidence_score_present = True
                        except ValueError:
                            confidence_score_present = False

                    if confidence_score_present:
                        # Handle rows with confidence score
                        confidence_score = float(cols[6])
                        if key not in target_positions:
                            target_positions[key] = line
                        else:
                            existing_line = target_positions[key]
                            existing_score = float(existing_line.split()[6])
                            if confidence_score > existing_score:
                                target_positions[key] = line
                    else:
                        # Handle rows without confidence score
                        if key not in target_positions:
                            target_positions[key] = line

                    # Categorize rows
                    if key in ref_positions:
                        matching_rows[key] = target_positions[key]
                        labeled_file.write(target_positions[key] + '\t1\n')
                    else:
                        non_matching_rows[key] = target_positions[key]
                        labeled_file.write(target_positions[key] + '\t0\n')
    except FileNotFoundError:
        print(f"Error: Target file '{target_file_path}' not found.")
        sys.exit(1)

    # Write results to output files
    matching_count = 0
    with open(matching_output_file, 'w') as match_file:
        for row in matching_rows.values():
            match_file.write(row + '\n')
            matching_count += 1

    non_matching_count = 0
    with open(non_matching_output_file, 'w') as non_match_file:
        for row in non_matching_rows.values():
            non_match_file.write(row + '\n')
            non_matching_count += 1

    print("reference file rows", len(ref_positions))
    print("target file rows", len(target_positions))
    ref_only_positions = set(ref_positions.keys()) - set(target_positions.keys())
    ref_only_rows = [ref_positions[pos] for pos in ref_only_positions]

    ref_only_count = 0
    with open(ref_only_output_file, 'w') as ref_only_file:
        for row in ref_only_rows:
            ref_only_file.write(row + '\n')
            ref_only_count += 1

    print(f"Matching rows written to: {matching_output_file} ({matching_count} rows)")
    print(f"Non-matching rows written to: {non_matching_output_file} ({non_matching_count} rows)")
    print(f"Reference-only rows written to: {ref_only_output_file} ({ref_only_count} rows)")
    print(f"File with labeled rows written to: {labeled_output_file}")

    # Calculate metrics
    TP = matching_count
    FP = non_matching_count
    FN = ref_only_count

    precision = TP / (TP + FP) if (TP + FP) != 0 else 0.0
    recall = TP / (TP + FN) if (TP + FN) != 0 else 0.0
    f1_score = (2 * precision * recall) / (precision + recall) if (precision + recall) != 0 else 0.0

    print(f"\nPrecision: {precision:.4f}")
    print(f"Recall:    {recall:.4f}")
    print(f"F1 Score:  {f1_score:.4f}")

if __name__ == "__main__":
    main()
