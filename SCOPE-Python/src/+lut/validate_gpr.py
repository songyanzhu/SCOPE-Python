import pandas as pd
import numpy as np
import os
import joblib # For loading GPR model
import warnings

# Mock external functions
from .plot_1to1 import plot_1to1

# Assumed external function to read Actot (simple wrapper for readtable)
def read_actot(path, var_name='Actot'):
    """Mock for lut.read_actot, reads the Actot column from a CSV."""
    df = pd.read_csv(path)
    if var_name not in df.columns:
        raise ValueError(f"Validation file {path} must contain column '{var_name}'.")
    return df[var_name].values


def validate_gpr(gpr_path=None, val_in_path=None, val_out_path=None):
    """
    // src/+lut/validate_gpr.m
    Validates a trained GPR model against a known validation set.
    """
    if gpr_path is None:
        in_dir = os.path.join('..', 'exercise')
        gpr_path = os.path.join(in_dir, 'gpr_Actot.joblib')
        val_in_path = os.path.join(in_dir, 'validation_in.csv')
        val_out_path = os.path.join(in_dir, 'validation_out.csv')
        
    fig_path = os.path.join(os.path.dirname(val_in_path), 'validation_gpr.png')
    
    assert os.path.exists(gpr_path), \
        f"Did not find `{gpr_path}` file.\nHave you trained the gaussian process regression (GPR) with lut.train_gpr()?"
    assert os.path.exists(val_in_path), \
        f"Did not find `{val_in_path}` file.\nHave you generated the validation input with lut.pick100()?"
    assert os.path.exists(val_out_path), \
        f"Did not find `{val_out_path}` file.\nHave you copied the `fluxes.csv` from `../output` after SCOPE run on validation_in.csv?"

    # --- 1. Load Model and Data ---
    try:
        model_obj = joblib.load(gpr_path)
        gprMdl = model_obj['gprMdl']
        scaler = model_obj['scaler']
        feature_cols = model_obj['features']
        var_name = model_obj['var_name']
    except Exception as e:
        raise RuntimeError(f"Failed to load GPR model and/or scaler: {e}")

    val_in = pd.read_csv(val_in_path)
    # MATLAB: val_out = lut.read_actot(val_out_path);
    val_out = read_actot(val_out_path, var_name)
    
    # --- 2. Predict ---
    start_time = time.time()
    
    X_val = val_in[feature_cols]
    X_val_scaled = scaler.transform(X_val)

    # MATLAB: res = predict(gprMdl, val_in);
    res, _ = gprMdl.predict(X_val_scaled, return_std=True)
    
    end_time = time.time()
    print(f"Prediction finished in {end_time - start_time:.2f} seconds.")
    
    # --- 3. Plot ---
    # MATLAB: lut.plot_1to1(val_out, res, 'GPR', fig_path)
    plot_1to1(val_out, res, 'GPR', fig_path, var_name)

if __name__ == '__main__':
    validate_gpr()