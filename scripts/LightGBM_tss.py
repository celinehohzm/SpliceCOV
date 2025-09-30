#!/usr/bin/env python3
import os, sys, argparse, pickle
import numpy as np
import pandas as pd
import lightgbm as lgb
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder

# -----------------------
# Path resolution (Option B)
# -----------------------
def _default_model_dir():
    # Highest priority: explicit env var
    mdir = os.environ.get("SPLICECOV_MODEL_DIR")
    if mdir:
        return os.path.abspath(mdir)
    # Next: installed helpers dir (launcher sets this)
    hdir = os.environ.get("SPLICECOV_HELPERS_DIR")
    if hdir:
        return os.path.abspath(os.path.join(hdir, "model_output"))
    # Fallback: next to this file
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "model_output"))

MODEL_DIR_DEFAULT = _default_model_dir()

# Default filenames (basenames)
TRAINING_FILE_BASENAME = "liver.processtiebrush_round2.tsstes.ptf_with_labels.txt"
TSS_MODEL_BASENAME     = "tss_lightgbm_model.txt"
CPAS_MODEL_BASENAME    = "cpas_lightgbm_model.txt"
TSS_ENC_BASENAME       = "tss_label_encoders.pkl"
CPAS_ENC_BASENAME      = "cpas_label_encoders.pkl"

# -----------------------
# I/O utils
# -----------------------
def _resolve(path_or_name, base):
    if not path_or_name:
        return None
    return path_or_name if os.path.isabs(path_or_name) else os.path.join(base, path_or_name)

def _safe_best_iter(model):
    # best_iteration_ is None when early stopping not triggered; LightGBM accepts None
    return getattr(model, "best_iteration", None)

def _save(obj, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump(obj, f)

def _load(path):
    with open(path, "rb") as f:
        return pickle.load(f)

# -----------------------
# Feature config
# -----------------------
FEATURES = ["value2", "value3", "value4", "value5", "value6", "value7", "value8"]
COLS_TRAIN = [
    "chromosome", "position", "unused", "row_type",
    "value1", "value2", "value3", "value4", "value5", "value6", "value7", "value8", "label"
]
COLS_TEST = [
    "chromosome", "position", "unused", "row_type",
    "value1", "value2", "value3", "value4", "value5", "value6", "value7", "value8"
]

def _prep_X(df, encoders):
    X = df[FEATURES].copy()
    # Normalize value3 by max from training
    vmax = encoders.get("value3_max", 1.0) or 1.0
    X["value3_norm"] = (X["value3"] / vmax).fillna(0.0)
    X.drop(columns=["value3"], inplace=True)
    # Encode categorical value8
    le = encoders["le_value8"]
    X["value8_encoded"] = le.transform(X["value8"].astype(str))
    X.drop(columns=["value8"], inplace=True)
    return X

def _fit_encoders(df):
    # compute encoders/statistics on TRAIN subset
    vmax = df["value3"].max()
    if pd.isna(vmax) or vmax == 0:
        vmax = 1.0
    le8 = LabelEncoder().fit(df["value8"].astype(str))
    return {"value3_max": float(vmax), "le_value8": le8}

# -----------------------
# Training
# -----------------------
def _train_one(df_full, row_type, model_path, enc_path, print_row_index=0):
    print(f"[train] row_type={row_type}")
    df = df_full[df_full["row_type"] == row_type].copy()
    if df.empty:
        print(f"[train] No rows for {row_type}; skipping.")
        return None

    # Basic class check
    pos = (df["label"] == 1).sum()
    neg = (df["label"] == 0).sum()
    if pos == 0 or neg == 0:
        print(f"[train] Not enough class diversity for {row_type} (pos={pos}, neg={neg}); skipping.")
        return None

    # Balance via undersampling (cap negatives at positives)
    if neg > pos:
        df_neg = df[df["label"] == 0].sample(n=pos, random_state=42)
        df_pos = df[df["label"] == 1]
        df = pd.concat([df_pos, df_neg]).sample(frac=1.0, random_state=42).reset_index(drop=True)

    # Encoders/statistics
    enc = _fit_encoders(df)
    _save(enc, enc_path)
    print(f"[train] saved encoders -> {enc_path}")

    # Features/labels
    X = _prep_X(df, enc)
    y = df["label"].astype(int)

    if 0 <= print_row_index < len(X):
        print(f"[train] sample features idx {print_row_index}: {X.iloc[print_row_index].to_dict()}")

    X_tr, X_te, y_tr, y_te = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )

    # LightGBM params
    params = {
        "objective": "binary",
        "metric": "binary_logloss",
        "verbosity": -1,
        "learning_rate": 0.05,
        "num_leaves": 31,
        "max_depth": -1,
        "min_data_in_leaf": 20,
    }

    dtr = lgb.Dataset(X_tr, label=y_tr)
    dte = lgb.Dataset(X_te, label=y_te)

    print(f"[train] training model -> {model_path}")
    model = lgb.train(
        params, dtr, valid_sets=[dte], num_boost_round=1000,
        callbacks=[lgb.early_stopping(stopping_rounds=50)]
    )
    model.save_model(model_path)
    print(f"[train] saved model -> {model_path}")
    return model

def train_models(training_path, tss_model, cpas_model, tss_enc, cpas_enc, print_row_index=0):
    print(f"[train] loading training data: {training_path}")
    df = pd.read_csv(training_path, sep="\t", header=None)
    if df.shape[1] != len(COLS_TRAIN):
        raise ValueError(f"Unexpected training cols: {df.shape[1]} (expected {len(COLS_TRAIN)})")
    df.columns = COLS_TRAIN

    _train_one(df, "TSS",  tss_model,  tss_enc,  print_row_index=print_row_index)
    _train_one(df, "CPAS", cpas_model, cpas_enc, print_row_index=print_row_index)

# -----------------------
# Scoring
# -----------------------
def _load_model_and_enc(model_path, enc_path):
    if not os.path.exists(model_path) or not os.path.exists(enc_path):
        return None, None
    model = lgb.Booster(model_file=model_path)
    enc = _load(enc_path)
    return model, enc

def score(testing_path, output_path, tss_model, cpas_model, tss_enc, cpas_enc, threshold=0.4):
    if not (0.0 <= float(threshold) <= 1.0):
        raise ValueError(f"threshold must be in [0,1], got {threshold}")
    print(f"[score] threshold={threshold}")

    df = pd.read_csv(testing_path, sep="\t", header=None)
    if df.shape[1] != len(COLS_TEST):
        raise ValueError(f"Unexpected testing cols: {df.shape[1]} (expected {len(COLS_TEST)})")
    df.columns = COLS_TEST

    # Prepare default outputs
    df["confidence_score"] = np.nan
    df["predicted_label"] = np.nan

    for rt, mpath, epath in (("TSS", tss_model, tss_enc), ("CPAS", cpas_model, cpas_enc)):
        print(f"[score] {rt}: loading model/encoders")
        model, enc = _load_model_and_enc(mpath, epath)
        if model is None or enc is None:
            print(f"[score] {rt}: model or encoders not found; skipping.")
            continue

        sub_idx = df.index[df["row_type"] == rt]
        if len(sub_idx) == 0:
            print(f"[score] {rt}: no rows; skipping.")
            continue

        X = _prep_X(df.loc[sub_idx], enc)
        y_prob = np.clip(model.predict(X, num_iteration=_safe_best_iter(model)), 0, 1)
        y_hat = (y_prob >= float(threshold)).astype(int)

        df.loc[sub_idx, "confidence_score"] = y_prob
        df.loc[sub_idx, "predicted_label"]  = y_hat

    # If any rows still NaN (no model), default to 0
    df["confidence_score"] = df["confidence_score"].fillna(0.0)
    df["predicted_label"]  = df["predicted_label"].fillna(0).astype(int)

    print(f"[score] writing -> {output_path}")
    df.to_csv(output_path, sep="\t", index=False)
    print("[score] done.")

# -----------------------
# CLI
# -----------------------
def main():
    p = argparse.ArgumentParser(description="Train and score LightGBM models for TSS/CPAS.")
    p.add_argument("-i", "--input",  dest="testing_file",  required=True,
                   help="Input file (TSSTES features without label).")
    p.add_argument("-o", "--output", dest="output_file",   required=True,
                   help="Output path for scored TSV.")
    p.add_argument("-s", "--threshold", type=float, default=0.4,
                   help="Prediction threshold in [0,1] (default: 0.4).")
    p.add_argument("--train", action="store_true",
                   help="Train models before scoring (or when models/encoders are missing).")
    p.add_argument("--model-dir", default=MODEL_DIR_DEFAULT,
                   help="Directory containing/holding models & encoders (default: resolved).")
    p.add_argument("--training-file", default=TRAINING_FILE_BASENAME,
                   help="Labeled training TSV (relative to model-dir unless absolute).")
    p.add_argument("--row", type=int, default=0, help="Row index to print during training.")
    args = p.parse_args()

    model_dir = os.path.abspath(args.model_dir)
    os.makedirs(model_dir, exist_ok=True)

    training_file = _resolve(args.training_file, model_dir)
    tss_model = _resolve(TSS_MODEL_BASENAME, model_dir)
    cpas_model = _resolve(CPAS_MODEL_BASENAME, model_dir)
    tss_enc   = _resolve(TSS_ENC_BASENAME,   model_dir)
    cpas_enc  = _resolve(CPAS_ENC_BASENAME,  model_dir)

    print(f"[paths] MODEL_DIR={model_dir}")
    print(f"[paths] TSS_MODEL={tss_model}")
    print(f"[paths] CPAS_MODEL={cpas_model}")
    print(f"[paths] TSS_ENC={tss_enc}")
    print(f"[paths] CPAS_ENC={cpas_enc}")

    need_train = args.train or (not (os.path.exists(tss_model) and os.path.exists(tss_enc))
                                or not (os.path.exists(cpas_model) and os.path.exists(cpas_enc)))
    if need_train:
        if not os.path.exists(training_file):
            print(f"[train] training requested/required but not found: {training_file}", file=sys.stderr)
            sys.exit(1)
        train_models(training_file, tss_model, cpas_model, tss_enc, cpas_enc, print_row_index=args.row)

    score(args.testing_file, args.output_file, tss_model, cpas_model, tss_enc, cpas_enc, threshold=args.threshold)

if __name__ == "__main__":
    main()
