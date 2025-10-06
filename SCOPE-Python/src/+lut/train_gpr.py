import pandas as pd
import numpy as np
import os
import warnings
import joblib # For saving/loading models (.mat/.joblib equivalent)

# Scikit-learn for GPR and cross-validation
try:
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
    from sklearn.model_selection import train_test_split
    from sklearn.preprocessing import StandardScaler
except ImportError:
    warnings.warn("Scikit-learn (needed for GPR) is not installed.")
    GaussianProcessRegressor = None

# Mock external functions
from .plot_1to1 import plot_1to1

def train_gpr(lut_in_path=None, lut_out_path=None, var_name='Actot'):
    """
    // src/+lut/train_gpr.m
    Trains a Gaussian Process Regression (GPR) model for parameter retrieval.
    """
    if GaussianProcessRegressor is None:
        raise ImportError("Scikit-learn is required for train_gpr.")

    if lut_in_path is None:
        in_dir = os.path.join('..', 'exercise')
        lut_in_path = os.path.join(in_dir, 'lut_in.csv')
        lut_out_path = os.path.join(in_dir, 'lut_out.csv')

    mat_out = os.path.join(os.path.dirname(lut_in_path), f'gpr_{var_name}.joblib')
    
    assert os.path.exists(lut_in_path), \
        f"Did not find `{lut_in_path}` file.\nHave you generated the LUT input with lut.generate_lut_input()?"
    assert os.path.exists(lut_out_path), \
        f"Did not find `{lut_out_path}` file.\nHave you copied the `fluxes.csv` from `../output` after SCOPE run on lut_in.csv?"
    
    lut_in = pd.read_csv(lut_in_path)
    
    # MATLAB: opt = detectImportOptions(lut_out_path); flu = readtable(lut_out_path, opt);
    flu = pd.read_csv(lut_out_path)
    
    # Align tables and extract input features (X) and target (y)
    if var_name not in flu.columns:
        raise ValueError(f"Variable '{var_name}' not found in LUT output file.")
        
    # Merge input parameters and SCOPE output
    # Assuming 't' is the time column used for merging/alignment if required, otherwise just appending is enough
    lut_data = pd.merge(lut_in, flu[['t', var_name]], on='t', how='inner')
    
    # Features (X): All input parameters except 't'
    X = lut_data.drop(columns=['t', var_name])
    y = lut_data[var_name]
    
    # --- Cross-Validation Split ---
    # MATLAB: cv = cvpartition(size(lut_in,1),'HoldOut',0.3); idx = cv.test;
    # MATLAB's HoldOut 0.3 means 30% test data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42
    )

    # --- Training (Standardize is handled explicitly in Scikit-learn) ---
    print("Started training gaussian process regression.\nUsually takes 1 minute.")
    
    # 1. Standardize the data (as requested by 'Standardize', 1)
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # 2. Define GPR Model (RBF kernel is common for environmental data)
    kernel = C(1.0) * RBF(1.0)
    gprMdl = GaussianProcessRegressor(
        kernel=kernel, 
        alpha=(y_train.std() / 10)**2, # Noise level estimation
        normalize_y=False, 
        random_state=42
    )
    
    # 3. Train
    gprMdl.fit(X_train_scaled, y_train.values)
    
    # --- Save Model and Scaler ---
    # MATLAB: save(mat_out, 'gprMdl') -> Python saves gprMdl and scaler
    model_obj = {'gprMdl': gprMdl, 'scaler': scaler, 'features': X.columns.tolist(), 'var_name': var_name}
    joblib.dump(model_obj, mat_out)
    print(f"GPR model is saved in `{mat_out}`")
    
    # --- Prediction and Validation Plot ---
    
    # MATLAB: res = predict(gprMdl, lut_test);
    res, std_dev = gprMdl.predict(X_test_scaled, return_std=True)
    
    fig_path = os.path.join(os.path.dirname(lut_in_path), f'{var_name}_gpr.png')
    # lut.plot_1to1(lut_test.(var_name), res, 'GPR', fig_path, var_name)
    plot_1to1(y_test.values, res, 'GPR', fig_path, var_name)