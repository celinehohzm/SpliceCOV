import pandas as pd
import lightgbm as lgb
import os
import pickle
import numpy as np
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split


TRAINING_FILE = './liver.processtiebrush_round2.tsstes.ptf_with_labels.txt'
TSS_MODEL_FILE = './model_output/tss_lightgbm_model.txt'
CPAS_MODEL_FILE = './model_output/cpas_lightgbm_model.txt'
TSS_ENCODERS_FILE = './model_output/tss_label_encoders.pkl'
CPAS_ENCODERS_FILE = './model_output/cpas_label_encoders.pkl'
OUTPUT_DIR = './model_output/'


# Function to train models
def train_model(model_type, print_row_index=0):
    print(f"Starting model training for {model_type}...")

    # Load the input file
    print(f"Loading data from '{TRAINING_FILE}'...")
    data = pd.read_csv(TRAINING_FILE, sep='\t', header=None)

    # Assign column names
    columns = ['chromosome', 'position', 'unused', 'row_type', 'value1', 'value2', 'value3',
               'value4', 'value5', 'value6', 'value7', 'value8', 'label']
    data.columns = columns

    # Filter data based on row type
    print(f"Filtering data for '{model_type}' rows...")
    data = data[data['row_type'] == model_type]

    # Separate positive and negative samples
    positive_samples = data[data['label'] == 1]
    negative_samples = data[data['label'] == 0]

    # Undersample negative samples to match the number of positive samples
    print("Balancing the dataset by undersampling negative samples...")
    negative_samples = negative_samples.sample(n=len(positive_samples), random_state=42)
    balanced_data = pd.concat([positive_samples, negative_samples]).sample(frac=1, random_state=42)

    # Check class distribution after balancing
    num_positive = balanced_data[balanced_data['label'] == 1].shape[0]
    num_negative = balanced_data[balanced_data['label'] == 0].shape[0]
    print(f"Number of positive samples: {num_positive}")
    print(f"Number of negative samples: {num_negative}")

    # Extract features and target variable
    print("Extracting features and target variable...")
    features = ['value2', 'value3', 'value4', 'value5', 'value6', 'value7', 'value8']
    X = balanced_data[features].copy()
    y = balanced_data['label']

    # Preprocess features
    print("Preprocessing features...")
    
    # Normalize numerical feature 'value3'
    print("Normalizing 'value3' to range [0, 1]...")
    value3_max = X['value3'].max()
    X['value3_norm'] = X['value3'] / value3_max
    X = X.drop('value3', axis=1)

    # Encode 'value8' using LabelEncoder
    print("Encoding 'value8'...")
    le_value8 = LabelEncoder()
    X['value8_encoded'] = le_value8.fit_transform(X['value8'])
    X = X.drop('value8', axis=1)

    # Save encoders and max values
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    encoders_file = TSS_ENCODERS_FILE if model_type == 'TSS' else CPAS_ENCODERS_FILE
    print(f"Saving encoders and max values to '{encoders_file}'...")
    encoders = {
        'le_value8': le_value8,
        'value3_max': value3_max
    }
    with open(encoders_file, 'wb') as f:
        pickle.dump(encoders, f)

    # Split the data
    print("\nSplitting data into training and testing sets...")
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )

    # Create LightGBM datasets
    print("Creating LightGBM datasets...")
    train_data = lgb.Dataset(X_train, label=y_train)
    test_data = lgb.Dataset(X_test, label=y_test)

    # Set up parameters
    print("Setting up LightGBM parameters...")
    params = {
        'objective': 'binary',
        'metric': 'binary_logloss',
        'verbosity': -1,
        'learning_rate': 0.05,
        'num_leaves': 31,
        'max_depth': -1,
        'min_data_in_leaf': 20
    }

    # Train the model with early stopping
    print(f"Training the LightGBM model for {model_type} rows...")
    model = lgb.train(
        params,
        train_data,
        valid_sets=[test_data],
        num_boost_round=1000,
        callbacks=[lgb.early_stopping(stopping_rounds=50)]
    )

    # Save the trained model
    model_file = TSS_MODEL_FILE if model_type == 'TSS' else CPAS_MODEL_FILE
    print(f"Saving the trained model to '{model_file}'...")
    model.save_model(model_file)

    print(f"{model_type} model training complete.")
    return model


# Function to score new data
def score_data(testing_file, output_file):
    print(f"Loading the testing data from '{testing_file}'...")
    data = pd.read_csv(testing_file, sep='\t', header=None)

    # Assign column names
    columns = ['chromosome', 'position', 'unused', 'row_type', 'value1', 'value2', 'value3',
               'value4', 'value5', 'value6', 'value7', 'value8']
    data.columns = columns

    predictions = []

    for row_type in ['TSS', 'CPAS']:
        print(f"Scoring {row_type} rows...")

        # Load the respective model and encoders
        model_file = TSS_MODEL_FILE if row_type == 'TSS' else CPAS_MODEL_FILE
        encoders_file = TSS_ENCODERS_FILE if row_type == 'TSS' else CPAS_ENCODERS_FILE
        model = lgb.Booster(model_file=model_file)
        with open(encoders_file, 'rb') as f:
            encoders = pickle.load(f)

        # Filter the rows for the current row_type
        subset = data[data['row_type'] == row_type]
        if subset.empty:
            continue

        # Extract features
        features = ['value2', 'value3', 'value4', 'value5', 'value6', 'value7', 'value8']
        X = subset[features].copy()

        # Preprocess features
        X['value3_norm'] = X['value3'] / encoders['value3_max']
        X = X.drop('value3', axis=1)
        X['value8_encoded'] = encoders['le_value8'].transform(X['value8'])
        X = X.drop('value8', axis=1)

        # Predict probabilities
        y_pred_prob = model.predict(X, num_iteration=model.best_iteration)

        # Assign predictions back to the original DataFrame
        data.loc[subset.index, 'confidence_score'] = y_pred_prob
        data.loc[subset.index, 'predicted_label'] = (y_pred_prob >= 0.4).astype(int)

    print(f"Saving scored data to '{output_file}'...")
    data.to_csv(output_file, sep='\t', index=False)


# Main execution
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Train and predict using TSS and CPAS models.')
    parser.add_argument('-i', '--input', dest='testing_file', required=True, help='Input file for prediction.')
    parser.add_argument('-o', '--output', dest='output_file', required=True, help='Output file for predictions.')

    args = parser.parse_args()

    # Train the models
    train_model('TSS')
    train_model('CPAS')

    # Score the data
    score_data(args.testing_file, args.output_file)
