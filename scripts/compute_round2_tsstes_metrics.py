#!/usr/bin/env python3

import sys
import pyBigWig
import numpy as np
import math

def main():
    if len(sys.argv) != 3:
        print("Usage: script.py <bundle_file> <bigwig_file>")
        sys.exit(1)
    bundle_file = sys.argv[1]
    bigwig_file = sys.argv[2]
    
    try:
        bw = pyBigWig.open(bigwig_file)
    except Exception as e:
        print(f"Error opening bigWig file: {e}")
        sys.exit(1)
    
    current_chr = None
    with open(bundle_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split()
            if cols[0] == 'bundle':
                current_chr = cols[1]
                print(line)
            elif cols[0] in ('tstart', 'tend'):
                position = int(cols[1])
                event_type = cols[0]
                # Now process this position
                # angle_metric, slope_metric = compute_metrics(bw, current_chr, position, event_type)
                angle_metric, slope_metric, smoothness_metric = compute_metrics(bw, current_chr, position, event_type)
                # Determine TRUE/FALSE based on slope and event type
                if (event_type == 'tstart' and slope_metric > 0) or (event_type == 'tend' and slope_metric < 0):
                    truth_value = "TRUE"
                else:
                    truth_value = "FALSE"
                # Output the original line with metrics and truth value appended
                # print(f"{line}\t{angle_metric:.2f}\t{slope_metric:.4f}\t{truth_value}")
                print(f"{line}\t{angle_metric:.2f}\t{slope_metric:.4f}\t{smoothness_metric:.4f}\t{truth_value}")
            else:
                # Other lines, just print as is
                print(line)
    bw.close()

# def compute_metrics(bw, chrom, position, event_type):
#     # Ensure the chromosome exists in the bigWig file
#     chrom_len = bw.chroms(chrom)
#     if chrom_len is None:
#         print(f"Chromosome {chrom} not found in bigWig file.")
#         sys.exit(1)
    
#     # Define the 50bp region for slope calculation
#     if event_type == 'tstart':
#         start = position
#         end = position + 50
#         if end > chrom_len:
#             end = chrom_len
#         positions = np.arange(start, end)
#     elif event_type == 'tend':
#         start = position - 50
#         if start < 0:
#             start = 0
#         end = position
#         positions = np.arange(start, end)
#     else:
#         positions = []
    
#     if len(positions) == 0:
#         slope_metric = 0.0
#         angle_metric = 0.0
#     else:
#         coverages = bw.values(chrom, positions[0], positions[-1] + 1)
#         coverages = np.nan_to_num(coverages)
#         x = positions - positions[0]
#         y = coverages
#         if len(x) < 2:
#             slope_metric = 0.0
#             angle_metric = 0.0
#         else:
#             # Fit a linear regression line to the coverage data
#             A = np.vstack([x, np.ones(len(x))]).T
#             m, c = np.linalg.lstsq(A, y, rcond=None)[0]
#             slope_metric = m
#             # Compute the angle of the slope in degrees
#             angle_rad = math.atan(m)
#             angle_metric = math.degrees(angle_rad)
#     return angle_metric, slope_metric

def compute_metrics(bw, chrom, position, event_type):
    # Ensure the chromosome exists in the bigWig file
    chrom_len = bw.chroms(chrom)
    if chrom_len is None:
        print(f"Chromosome {chrom} not found in bigWig file.")
        sys.exit(1)
    
    # Define the 50bp region for slope calculation
    if event_type == 'tstart':
        start = position
        end = position + 50
        if end > chrom_len:
            end = chrom_len
        positions = np.arange(start, end)
    elif event_type == 'tend':
        start = position - 50
        if start < 0:
            start = 0
        end = position
        positions = np.arange(start, end)
    else:
        positions = []
    
    if len(positions) == 0:
        slope_metric = 0.0
        angle_metric = 0.0
        smoothness_metric = 0.0
    else:
        coverages = bw.values(chrom, positions[0], positions[-1] + 1)
        coverages = np.nan_to_num(coverages)
        x = positions - positions[0]
        y = coverages
        
        if len(x) < 2:
            slope_metric = 0.0
            angle_metric = 0.0
            smoothness_metric = 0.0
        else:
            # Fit a linear regression line to the coverage data
            A = np.vstack([x, np.ones(len(x))]).T
            m, c = np.linalg.lstsq(A, y, rcond=None)[0]
            slope_metric = m
            angle_rad = math.atan(m)
            angle_metric = math.degrees(angle_rad)
            
            # Compute smoothness metric
            diffs = np.diff(y)
            mean_coverage = np.mean(y) + 1e-6  # Avoid division by zero
            relative_diffs = diffs / mean_coverage  # Normalize differences by mean coverage
            
            # Relative metrics
            mean_abs_relative_diff = np.mean(np.abs(relative_diffs)) + 1e-6
            std_relative_diffs = np.std(relative_diffs)
            cv_relative = std_relative_diffs / mean_abs_relative_diff  # Relative CV
            max_relative_diff = np.max(np.abs(relative_diffs))
            
            # Weighting factors
            alpha = 0.5  # CV weight
            beta = 1.0   # Max relative diff weight
            
            # Final smoothness metric
            smoothness_metric = 1 / (1 + alpha * cv_relative + beta * max_relative_diff)
            
            # # Debug output for specific position
            # if chrom == "chr1" and position == 13989:
            #     print(f"--- Debugging Smoothness Metric for {chrom}:{position} ---")
            #     print(f"Coverages: {y}")
            #     print(f"Differences: {diffs}")
            #     print(f"Mean Coverage: {mean_coverage}")
            #     print(f"Relative Differences: {relative_diffs}")
            #     print(f"Mean Absolute Relative Difference: {mean_abs_relative_diff}")
            #     print(f"Standard Deviation of Relative Differences: {std_relative_diffs}")
            #     print(f"Coefficient of Variation (CV Relative): {cv_relative}")
            #     print(f"Maximum Relative Difference: {max_relative_diff}")
            #     print(f"Smoothness Metric: {smoothness_metric}")
    
    return angle_metric, slope_metric, smoothness_metric






if __name__ == '__main__':
    main()
