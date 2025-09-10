#!/usr/bin/env python3

import pandas as pd
import lightgbm as lgb
import os
import pickle
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, accuracy_score, classification_report, confusion_matrix
import argparse
import sys

# Paths to the saved model and encoders
TRAINING_FILE = '/ccb/salz1/choh1/spliceCov/scripts/model_output/spleen.process_tiebrush.processed_round1.bundles_with_labels.txt'
MODEL_FILE    = '/ccb/salz1/choh1/spliceCov/scripts/model_output/0605_spleen_lightgbm_model.txt'
ENCODERS_FILE = '/ccb/salz1/choh1/spliceCov/scripts/model_output/0605_spleen_label_encoders.pkl'
FEATURE_IMPORTANCES_FILE = '/ccb/salz1/choh1/spliceCov/scripts/model_output/0605_spleen_norm_feature_importances.csv'
OUTPUT_DIR    = '/ccb/salz1/choh1/spliceCov/scripts/model_output/'

# Proximity filter: keep only one row per cluster within `window` bp on same chr, based on max score_col
def filter_by_proximity(df, window=5, score_col='num_samples'):
    # Sort by chromosome and position
    df_sorted = df.sort_values(['chromosome', 'position']).reset_index(drop=True)
    keep_indices = []
    prev_chr = None
    prev_pos = None
    max_idx  = None
    max_score = -np.inf
    for idx, row in df_sorted.iterrows():
        chrom = row['chromosome']
        pos   = row['position']
        score = row[score_col]
        if prev_chr == chrom and (pos - prev_pos) <= window:
            # same cluster
            if score > max_score:
                max_score = score
                max_idx   = idx
        else:
            # new cluster: save previous
            if max_idx is not None:
                keep_indices.append(max_idx)
            # start new cluster
            prev_chr   = chrom
            prev_pos   = pos
            max_score  = score
            max_idx    = idx
    # append last
    if max_idx is not None:
        keep_indices.append(max_idx)
    return df_sorted.loc[keep_indices].reset_index(drop=True)


def train_model(print_row_index=0, scaling_factor=1.0):
    # ... (unchanged training code remains here) ...
    pass  # omitted for brevity


def score_data(testing_file, output_file):
    # Load model and encoders
    model = lgb.Booster(model_file=MODEL_FILE)
    with open(ENCODERS_FILE, 'rb') as f:
        encoders = pickle.load(f)

    # Load data
    data = pd.read_csv(testing_file, sep='\t', header=None)
    # assign columns...
    if data.shape[1] == 13:
        data.columns = ['chromosome','position','junction_id','num_samples','strand',
                        'perc','cov_diff','perc_cov_diff','junc_len','smooth_metric',
                        'cov_change_dir','event','label']
    else:
        data.columns = ['chromosome','position','junction_id','num_samples','strand',
                        'perc','cov_diff','perc_cov_diff','junc_len','smooth_metric',
                        'cov_change_dir','event']

    # preserve raw
    smooth_metric_raw = data['smooth_metric'].copy()

    # build features
    X = data[['num_samples','perc','cov_diff','perc_cov_diff','junc_len','smooth_metric','cov_change_dir']].copy()
    X['perc_binarized'] = X['perc'].apply(lambda x: 1 if x == "1.0000-1.0000-1.0000-1.0000" else 0)
    X = X.drop('perc', axis=1)
    # num_samples norm
    cov_med = encoders.get('coverage_median',1.0)
    sf     = encoders.get('scaling_factor',1.0)
    X['num_samples_scaled'] = X['num_samples']/cov_med/sf
    X = X.drop('num_samples', axis=1)
    # keep raw junc_len and smooth_metric
    X['junc_len_log']      = np.log1p(X['junc_len'])  # you can remove log1p if raw needed
    X = X.drop('junc_len', axis=1)
    X['smooth_metric_scaled'] = (X['smooth_metric'] - encoders.get('smooth_metric_mean',0.0))/encoders.get('smooth_metric_std',1.0)
    X = X.drop('smooth_metric', axis=1)

    # predict
    y_pred_prob = model.predict(X, num_iteration=model.best_iteration)
    y_pred_prob = np.clip(y_pred_prob, 0, 1)
    y_pred      = (y_pred_prob >= 0.4).astype(int)

    # special override
    smooth_override_mask = (
        (smooth_metric_raw < 10) &
        (data['perc_cov_diff'] <= 0.25) &
        (data['junc_len'] <= 15)
    )
    y_pred_prob[smooth_override_mask] = 0.0
    y_pred[smooth_override_mask]      = 0
    # cov_change_dir override
    if 'cov_change_dir' in data.columns:
        mask = data['cov_change_dir']==0
        y_pred_prob[mask] = 0.0
        y_pred[mask]      = 0

    data['confidence_score'] = y_pred_prob
    data['predicted_label']  = y_pred

    # Apply proximity filter on scored output using num_samples (col 4)
    # ensure 'num_samples' exists in data
    df_filtered = filter_by_proximity(data, window=5, score_col='num_samples')

    # Save filtered
    df_filtered.to_csv(output_file, sep='\t', index=False)
    print(f"Filtered scored data saved to '{output_file}'.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input',dest='input_file', required=True)
    parser.add_argument('-o','--output',dest='output_file', required=True)
    parser.add_argument('-r','--row',dest='row_index',type=int, default=0)
    parser.add_argument('-s','--scale',dest='scaling_factor',type=float, default=1.0)
    args = parser.parse_args()

    if not os.path.exists(MODEL_FILE) or not os.path.exists(ENCODERS_FILE):
        train_model(print_row_index=args.row_index, scaling_factor=args.scaling_factor)
    if os.path.abspath(args.input_file) == os.path.abspath(TRAINING_FILE):
        train_model(print_row_index=args.row_index, scaling_factor=args.scaling_factor)
    score_data(args.input_file, args.output_file)
