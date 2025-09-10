#!/usr/bin/env python3
import sys
import pyBigWig
import numpy as np

SMALL_DELTA = 5   # for left/right mean windows
W = 50            # for smoothness window (±W around pos)

def window_vals(bw, chrom, pos, chrom_len, w=W):
    """Return coverage array for [pos-w, pos+w] ∩ [0, chrom_len), and index of pos within it."""
    start = max(0, pos - w)
    end   = min(chrom_len, pos + w + 1)  # end-exclusive
    arr = np.array(bw.values(chrom, start, end), dtype=np.float32)
    np.nan_to_num(arr, copy=False)
    center = pos - start  # index of pos inside arr
    return arr, center, start, end

def left_right_means(arr, center, strand, d=SMALL_DELTA):
    """Compute mean coverage on left/right d bases relative to center, respecting strand flip."""
    # raw left/right relative to genomic coordinate
    left  = arr[max(0, center - d): center]
    right = arr[center + 1: center + 1 + d]
    left_mean  = float(left.mean())  if left.size  else 0.0
    right_mean = float(right.mean()) if right.size else 0.0

    # keep the same semantics as your original function:
    # for '-' we conceptually flip sides when interpreting "left"/"right"
    if strand == '-':
        left_mean, right_mean = right_mean, left_mean

    # “coverage at pos” used only conceptually; compute in-case you want it later
    cov_at_pos = float(arr[center]) if 0 <= center < arr.size else 0.0

    # percent change and absolute change (same formulas you used)
    if left_mean > 0 and right_mean > 0:
        ratio = right_mean / left_mean if right_mean >= left_mean else left_mean / right_mean
        perc_change = 1.0 - min(ratio, 1.0/ratio) if ratio != 0 else 0.0
    else:
        perc_change = 0.0
    abs_change = abs(right_mean - left_mean)

    return cov_at_pos, left_mean, right_mean, perc_change, abs_change

def smoothness_from_window(arr, center, which):
    """
    Compute smoothness using diffs from a directional window:
      - JSTART: use coverage[pos : pos+W]
      - JEND  : use coverage[pos-W : pos]
    Exclude the junction’s own diff (first for JSTART, last for JEND).
    Returns (second_largest_abs_diff).
    """
    if which == 'JSTART':
        sub = arr[center : min(arr.size, center + W + 1)]
        diffs = np.diff(sub)
        others = diffs[1:] if diffs.size > 1 else np.empty(0, dtype=np.float32)
    else:  # 'JEND'
        sub = arr[max(0, center - W) : center + 1]
        diffs = np.diff(sub)
        others = diffs[:-1] if diffs.size > 1 else np.empty(0, dtype=np.float32)

    if others.size < 2:
        return 0.0
    a = np.abs(others)
    # second largest via partition (O(n))
    idx = -2
    part = np.partition(a, idx)
    return float(part[idx])

def process_junctions(bw_file, junc_file, small_delta=SMALL_DELTA, w=W):
    try:
        bw = pyBigWig.open(bw_file)
    except Exception as e:
        print(f"Error opening BigWig file '{bw_file}': {e}", file=sys.stderr)
        sys.exit(1)

    # cache lengths once
    try:
        chrom_len_map = bw.chroms()
    except Exception as e:
        print(f"Error reading chrom sizes: {e}", file=sys.stderr)
        sys.exit(1)

    try:
        with open(junc_file, 'r') as jf:
            for line in jf:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 7:
                    # malformed line (keep your warning behavior)
                    print(f"Warning: Skipping malformed line: {line.strip()}", file=sys.stderr)
                    continue

                chrom, start, end, name, cov, strand, percs = fields[:7]
                try:
                    start = int(start); end = int(end); cov = int(cov)
                except ValueError:
                    print(f"Warning: Non-integer values found in line: {line.strip()}", file=sys.stderr)
                    continue

                chrom_len = chrom_len_map.get(chrom)
                if chrom_len is None:
                    print(f"Chromosome {chrom} not found in BigWig file.", file=sys.stderr)
                    continue

                jlen = end - start

                # ------- START side -------
                arrS, cS, _, _ = window_vals(bw, chrom, start, chrom_len, w)
                cov_pos_S, lS, rS, pcS, abschgS = left_right_means(arrS, cS, strand, small_delta)
                secS = smoothness_from_window(arrS, cS, 'JSTART')
                # NOTE: your original 'smoothness_metric = abs(change_at_pos)/second_largest or abs(change_at_pos)'
                # here change_at_pos ~ abschgS (magnitude difference across sides near junction)
                smS = abschgS if secS == 0.0 else abs(abschgS / secS)

                # your original directional flag logic:
                jstart_flag = "TRUE" if ((strand == "+" and (rS - lS) < 0) or
                                          (strand == "-" and (rS - lS) > 0)) else "FALSE"
                encS = 1 if jstart_flag == "TRUE" else 0

                sys.stdout.write(
                    f"{chrom}\t{start}\t{name}\t{cov}\t{strand}\t{percs}\t{jlen}"
                    f"\t{pcS:.4f}\t{abs(abschgS):.0f}\t{smS:.4f}\t{encS}\tJSTART\n"
                )

                # ------- END side -------
                arrE, cE, _, _ = window_vals(bw, chrom, end, chrom_len, w)
                cov_pos_E, lE, rE, pcE, abschgE = left_right_means(arrE, cE, strand, small_delta)
                secE = smoothness_from_window(arrE, cE, 'JEND')
                smE = abschgE if secE == 0.0 else abs(abschgE / secE)

                jend_flag = "TRUE" if ((strand == "+" and (rE - lE) > 0) or
                                        (strand == "-" and (rE - lE) < 0)) else "FALSE"
                encE = 1 if jend_flag == "TRUE" else 0

                sys.stdout.write(
                    f"{chrom}\t{end+1}\t{name}\t{cov}\t{strand}\t{percs}\t{jlen}"
                    f"\t{pcE:.4f}\t{abs(abschgE):.0f}\t{smE:.4f}\t{encE}\tJEND\n"
                )

    except FileNotFoundError:
        print(f"Error: Junctions file '{junc_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing junctions file: {e}", file=sys.stderr)
        sys.exit(1)
    finally:
        try:
            bw.close()
        except Exception:
            pass

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python process_tiebrush_round1_juncs_splicecov.py <bigwig_file> <junctions_file>")
        sys.exit(1)
    process_junctions(sys.argv[1], sys.argv[2])
