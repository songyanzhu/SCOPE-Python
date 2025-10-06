import pandas as pd
import numpy as np
import os
import warnings

# Mock external functions
from .lut_search import lut_search
from .plot_1to1 import plot_1to1
from .validate_gpr import read_actot # Reuse helper from validate_gpr

def validate_rmse(lut_in_path=None, lut_out_path=None, val_in_path=None, val_out_path=None):
    """
    // src/+lut/validate_rmse.m
    Validates the LUT search (RMSE matching) method against a known validation set.
    """
    var_name = 'Actot'
    
    if lut_in_path is None:
        in_dir = os.path.join('..', 'exercise')
        lut_in_path = os.path.join(in_dir, 'lut_in.csv')
        lut_out_path = os.path.join(in_dir, 'lut_out.csv')
        val_in_path = os.path.join(in_dir, 'validation_in.csv')
        val_out_path = os.path.join(in_dir, 'validation_out.csv')
        
    fig_path = os.path.join(os.path.dirname(val_in_path), 'validation_rmse.png')
    
    assert os.path.exists(lut_in_path), \
        f"Did not find `{lut_in_path}` file.\nHave you generated the LUT input with lut.generate_lut_input()?"
    assert os.path.exists(lut_out_path), \
        f"Did not find `{lut_out_path}` file.\nHave you copied the `fluxes.csv` from `../output` after SCOPE run on lut_in.csv?"
    assert os.path.exists(val_in_path), \
        f"Did not find `{val_in_path}` file.\nHave you generated the validation input with lut.pick100()?"
    assert os.path.exists(val_out_path), \
        f"Did not find `{val_out_path}` file.\nHave you copied the `fluxes.csv` from `../output` after SCOPE run on validation_in.csv?"

    # --- 1. Load Data ---
    lut_in = pd.read_csv(lut_in_path)
    
    # MATLAB: lut_out = lut.read_actot(lut_out_path);
    lut_out = read_actot(lut_out_path, var_name)
    
    val_in = pd.read_csv(val_in_path)
    # MATLAB: val_out = lut.read_actot(val_out_path);
    val_out = read_actot(val_out_path, var_name)
    
    # --- 2. LUT Search ---
    # [res, res_std] = lut.lut_search(val_in, lut_in, lut_out);
    res, res_std = lut_search(val_in, lut_in, lut_out)
    
    # --- 3. Plot ---
    # lut.plot_1to1(val_out, res, 'RMSE', fig_path)
    plot_1to1(val_out, res, 'RMSE', fig_path, var_name)

if __name__ == '__main__':
    validate_rmse()