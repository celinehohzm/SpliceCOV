#!/usr/bin/env python3
import sys
import pandas as pd

def print_percentiles(data, col_name):
    """
    Given a pandas Series `data`, compute and print the 5th, 10th, ..., 95th percentiles.
    """
    # Define the desired quantile probabilities
    probs = [i/100 for i in range(5, 100, 5)]  # 0.05, 0.10, ..., 0.95
    percentiles = data.quantile(probs)
    print(f"Column {col_name} percentiles:")
    for p, val in zip(range(5, 100, 5), percentiles):
        print(f"  {p}th percentile: {val}")
    print()

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <sorted_input_file>", file=sys.stderr)
        sys.exit(1)

    input_file = sys.argv[1]

    # Read the file as a tab-delimited file with no header
    try:
        df = pd.read_csv(input_file, sep='\t', header=None)
    except Exception as e:
        print(f"Error reading '{input_file}': {e}", file=sys.stderr)
        sys.exit(1)

    # Columns of interest (1-based indices): 4, 7, 8, 9, 10, 13
    # Convert to 0-based indices for pandas:
    cols_to_analyze = {
        4: 3,
        7: 6,
        8: 7,
        9: 8,
        10: 9,
        13: 12
    }

    for col_num, idx in cols_to_analyze.items():
        if idx >= df.shape[1]:
            print(f"Error: Column {col_num} (index {idx}) not found in '{input_file}'", file=sys.stderr)
            continue

        # Ensure data is numeric
        try:
            series = pd.to_numeric(df.iloc[:, idx], errors='coerce').dropna()
        except Exception as e:
            print(f"Error converting column {col_num} to numeric: {e}", file=sys.stderr)
            continue

        print_percentiles(series, col_num)

if __name__ == "__main__":
    main()
