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
TRAINING_FILE = 'model_output/spleen.process_tiebrush.processed_round1.bundles_with_labels.txt'
MODEL_FILE = 'model_output/0606_spleen_no_normscale_lightgbm_model.txt'
ENCODERS_FILE = 'model_output/0606_spleen_no_normscale_label_encoders.pkl'
FEATURE_IMPORTANCES_FILE = 'model_output/0606_spleen_no_normscale_norm_feature_importances.csv'
OUTPUT_DIR = './model_output/'

def train_model(print_row_index=0, scaling_factor=1.0):
    print("Starting model training...")

    # Load the input file
    print(f"Loading data from '{TRAINING_FILE}'...")
    try:
        data = pd.read_csv(TRAINING_FILE, sep='\t', header=None)
    except Exception as e:
        print(f"Error loading training file: {e}", file=sys.stderr)
        sys.exit(1)

    # Define column names for the training file (13 columns: 12 features + label)
    columns_training = [
        'chromosome', 'position', 'junction_id', 'num_samples', 'strand',
        'perc', 'cov_diff', 'perc_cov_diff', 'junc_len', 'smooth_metric',
        'cov_change_dir', 'event', 'label'
    ]
    data.columns = columns_training

    # Check if 'label' column exists
    if 'label' not in data.columns:
        print("Error: 'label' column is missing in the training data.", file=sys.stderr)
        sys.exit(1)

    # Check class distribution
    num_positive = data[data['label'] == 1].shape[0]
    num_negative = data[data['label'] == 0].shape[0]
    print(f"Number of positive samples: {num_positive}")
    print(f"Number of negative samples: {num_negative}")

    if num_positive == 0:
        print("Error: No positive samples found in the dataset.", file=sys.stderr)
        sys.exit(1)

    # Calculate scale_pos_weight for handling class imbalance
    scale_pos_weight = num_negative / num_positive
    print(f"Scale_pos_weight: {scale_pos_weight}")

    # Extract features and target
    print("Extracting features and target variable...")
    X = data[['num_samples', 'perc', 'cov_diff', 'perc_cov_diff', 'junc_len', 'smooth_metric', 'cov_change_dir']].copy()
    y = data['label']

    # 1) Binarize 'perc'
    print("Binarizing 'perc' into two categories...")
    X['perc_binarized'] = X['perc'].apply(
        lambda x: 1 if x == "1.0000-1.0000-1.0000-1.0000" else 0
    )
    X = X.drop('perc', axis=1)

    # 2) Normalize 'num_samples' by median, then scale it (unchanged)
    print("Normalizing 'num_samples' by median, then scaling it...")
    coverage_median = X['num_samples'].median()
    if coverage_median == 0:
        print("Warning: median coverage is 0, skipping median normalization.")
        coverage_median = 1.0
    X['num_samples_med_norm'] = X['num_samples'] / coverage_median
    X['num_samples_scaled'] = X['num_samples_med_norm'] / scaling_factor
    X = X.drop(['num_samples', 'num_samples_med_norm'], axis=1)

    # 3) Keep 'junc_len' raw (no log-transform)
    X['junc_len_raw'] = X['junc_len']
    X = X.drop('junc_len', axis=1)

    # 4) Keep 'smooth_metric' raw (no standard scaling)
    X['smooth_metric_raw'] = X['smooth_metric']
    X = X.drop('smooth_metric', axis=1)

    print(f"Training features: {list(X.columns)}")

    # Save encoders/normalization parameters for future use
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    print(f"Saving encoders and normalization parameters to '{ENCODERS_FILE}'...")
    encoders = {
        'coverage_median': coverage_median,
        'scaling_factor': scaling_factor
    }
    with open(ENCODERS_FILE, 'wb') as f:
        pickle.dump(encoders, f)

    print(f"\n--- Features for Training Data Row Index: {print_row_index} ---")
    if 0 <= print_row_index < len(X):
        sample_features = X.iloc[print_row_index]
        print(sample_features.to_dict())
    else:
        print(f"Requested row index {print_row_index} is out of bounds. Total rows: {len(X)}")

    print("\nSplitting data into training and testing sets...")
    try:
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42, stratify=y
        )
    except Exception as e:
        print(f"Error during train-test split: {e}", file=sys.stderr)
        sys.exit(1)

    print("Creating LightGBM datasets...")
    train_data = lgb.Dataset(X_train, label=y_train)
    test_data = lgb.Dataset(X_test, label=y_test)

    print("Setting up LightGBM parameters...")
    params = {
        'objective': 'binary',
        'metric': 'binary_logloss',
        'scale_pos_weight': scale_pos_weight,
        'verbosity': -1,
        'learning_rate': 0.05,
        'num_leaves': 31,
        'max_depth': -1,
        'max_bin': 512,
        'min_data_in_leaf': 10
    }

    print("Training the LightGBM model with early stopping...")
    try:
        model = lgb.train(
            params,
            train_data,
            valid_sets=[test_data],
            num_boost_round=1000,
            callbacks=[lgb.early_stopping(stopping_rounds=50)]
        )
    except Exception as e:
        print(f"Error during model training: {e}", file=sys.stderr)
        sys.exit(1)

    print("Evaluating the model on the test set...")
    try:
        y_pred_prob = model.predict(X_test, num_iteration=model.best_iteration)
        threshold = 0.5
        y_pred = (y_pred_prob >= threshold).astype(int)
        auc = roc_auc_score(y_test, y_pred_prob)
        accuracy = accuracy_score(y_test, y_pred)
        print(f'Model training complete. AUC: {auc:.4f}, Accuracy: {accuracy:.4f}')
    except Exception as e:
        print(f"Error during model evaluation: {e}", file=sys.stderr)
        sys.exit(1)

    print("\nClassification Report:")
    print(classification_report(y_test, y_pred))
    print("Confusion Matrix:")
    print(confusion_matrix(y_test, y_pred))

    print("\nFeature importances:")
    feature_importances = pd.DataFrame({
        'feature': model.feature_name(),
        'importance': model.feature_importance(importance_type='gain')
    })
    feature_importances_sorted = feature_importances.sort_values(by='importance', ascending=False)
    for _, row in feature_importances_sorted.iterrows():
        print(f"  {row['feature']}: {row['importance']:.2f}")

    print(f"\nSaving feature importances to '{FEATURE_IMPORTANCES_FILE}'...")
    try:
        feature_importances_sorted.to_csv(FEATURE_IMPORTANCES_FILE, index=False)
        print("Feature importances saved successfully.")
    except Exception as e:
        print(f"Error saving feature importances: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"\nSaving the trained model to '{MODEL_FILE}'...")
    try:
        model.save_model(MODEL_FILE)
    except Exception as e:
        print(f"Error saving the model: {e}", file=sys.stderr)
        sys.exit(1)

    print("Training complete.")
    return model, encoders


def score_data(testing_file, output_file, threshold=0.4):
    # --- NEW: threshold is now a parameter (default 0.4) ---
    if not (0.0 <= float(threshold) <= 1.0):
        print(f"Error: threshold must be in [0,1], got {threshold}", file=sys.stderr)
        sys.exit(1)
    print(f"Using classification threshold: {threshold}")

    print(f"Loading the pretrained model from '{MODEL_FILE}'...")
    try:
        model = lgb.Booster(model_file=MODEL_FILE)
    except Exception as e:
        print(f"Error loading the model: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading encoders and normalization parameters from '{ENCODERS_FILE}'...")
    try:
        with open(ENCODERS_FILE, 'rb') as f:
            encoders = pickle.load(f)
    except Exception as e:
        print(f"Error loading encoders: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading input data from '{testing_file}'...")
    try:
        data = pd.read_csv(testing_file, sep='\t', header=None)
    except Exception as e:
        print(f"Error loading testing file: {e}", file=sys.stderr)
        sys.exit(1)

    columns_training = [
        'chromosome','position','junction_id','num_samples','strand','perc',
        'cov_diff','perc_cov_diff','junc_len','smooth_metric','cov_change_dir','event','label'
    ]
    columns_testing = [
        'chromosome','position','junction_id','num_samples','strand','perc',
        'cov_diff','perc_cov_diff','junc_len','smooth_metric','cov_change_dir','event'
    ]

    if data.shape[1] == 13:
        data.columns = columns_training
    elif data.shape[1] == 12:
        data.columns = columns_testing
    else:
        print(f"Error: Unexpected number of columns: {data.shape[1]}", file=sys.stderr)
        sys.exit(1)

    print("Extracting features...")
    smooth_metric_raw = data['smooth_metric'].copy()
    X = data[['num_samples','perc','cov_diff','perc_cov_diff','junc_len','smooth_metric','cov_change_dir']].copy()

    print("Binarizing 'perc' into two categories...")
    X['perc_binarized'] = X['perc'].apply(
        lambda x: 1 if x == "1.0000-1.0000-1.0000-1.0000" else 0
    )
    X = X.drop('perc', axis=1)

    coverage_median = encoders.get('coverage_median', 1.0)
    scaling_factor = encoders.get('scaling_factor', 1.0)
    if coverage_median == 0:
        coverage_median = 1.0
    print(f"Normalizing 'num_samples' by median={coverage_median}, then scaling by dividing by {scaling_factor}...")
    X['num_samples_med_norm'] = X['num_samples'] / coverage_median
    X['num_samples_scaled'] = X['num_samples_med_norm'] / scaling_factor
    X = X.drop(['num_samples','num_samples_med_norm'], axis=1)

    X['junc_len_raw'] = X['junc_len']; X = X.drop('junc_len', axis=1)
    X['smooth_metric_raw'] = X['smooth_metric']; X = X.drop('smooth_metric', axis=1)

    print(f"Prediction features: {list(X.columns)}")

    print("Making predictions...")
    try:
        y_pred_prob = model.predict(X, num_iteration=model.best_iteration)
    except Exception as e:
        print(f"Error during prediction: {e}", file=sys.stderr)
        sys.exit(1)

    y_pred_prob = np.clip(y_pred_prob, 0, 1)

    # --- Threshold now configurable ---
    y_pred = (y_pred_prob >= float(threshold)).astype(int)

    # =========================== SPECIAL OVERRIDE ===========================
    smooth_override_mask = (smooth_metric_raw < 15) & (data['perc_cov_diff'] <= 0.25) & (data['junc_len'] <= 20)
    y_pred_prob[smooth_override_mask] = 0.0
    y_pred[smooth_override_mask] = 0
    # ======================================================================

    if 'cov_change_dir' in data.columns:
        override_mask_cov_change = (data['cov_change_dir'] == 0)
        y_pred_prob[override_mask_cov_change] = 0.0
        y_pred[override_mask_cov_change] = 0

    data['confidence_score'] = y_pred_prob
    data['predicted_label'] = y_pred

    print(f"Saving the results to '{output_file}'...")
    try:
        data.to_csv(output_file, sep='\t', index=False)
        print(f"Scored data saved to '{output_file}'.")
    except Exception as e:
        print(f"Error saving scored data: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Train a LightGBM model and score new data.')
    parser.add_argument('-i', '--input', dest='input_file', required=True,
                        help='Path to the input file (training or testing).')
    parser.add_argument('-o', '--output', dest='output_file', required=True,
                        help='Path to save the output file (scores).')
    parser.add_argument('-r', '--row', dest='row_index', type=int, default=0,
                        help='Row index to print features for during training (default: 0).')
    parser.add_argument('-x', '--scale', dest='scaling_factor', type=float, default=1.0,
                        help='Scaling factor to normalize num_samples by dividing by this number (default: 1.0).')
    # --- NEW: threshold flag ---
    parser.add_argument('-s', '--threshold', dest='threshold', type=float, default=0.4,
                        help='Prediction threshold in [0,1] for converting probabilities to labels (default: 0.4).')

    args = parser.parse_args()

    # Check if model and encoders exist
    if not os.path.exists(MODEL_FILE) or not os.path.exists(ENCODERS_FILE):
        print("Pretrained model or encoders not found. Initiating training...")
        train_model(print_row_index=args.row_index, scaling_factor=args.scaling_factor)
    else:
        print("Pretrained model and encoders found. Skipping training.")

    # If the input file is the same as the training file, retrain the model.
    if os.path.abspath(args.input_file) == os.path.abspath(TRAINING_FILE):
        print("Input file is the training file. Retraining the model...")
        train_model(print_row_index=args.row_index, scaling_factor=args.scaling_factor)

    # Score data with user-specified threshold
    score_data(args.input_file, args.output_file, threshold=args.threshold)
    print("Done.")
