import pandas as pd
import numpy as np
import os
import warnings
import joblib # For loading GPR model

# Mock external functions
from .csv2image_plane import csv2image_plane
from .plot_image import plot_image 

# Assumed external function to write image without georeference
def write_tiff(im, out_path):
    """Mock for lut.write_tiff - saves a numpy array as a GeoTIFF (simplified)."""
    try:
        import rasterio
        profile = {'driver': 'GTiff', 'dtype': 'float32', 'nodata': np.nan, 
                   'width': im.shape[1], 'height': im.shape[0], 'count': 1,
                   'compress': 'PACKBITS'}
        with rasterio.open(out_path, 'w', **profile) as dst:
            dst.write(im, 1)
    except ImportError:
        warnings.warn("Rasterio not found. Cannot write GeoTIFF.")

def use_gpr(gpr_path=None, full_set_path=None):
    """
    // src/+lut/use_gpr.m
    Uses the trained GPR model to predict a response variable on a full dataset.
    """
    if gpr_path is None:
        in_dir = os.path.join('..', 'exercise')
        gpr_path = os.path.join(in_dir, 'gpr_Actot.joblib') # Assuming a default name
        full_set_path = os.path.join(in_dir, 'full_set.csv')

    out_dir = os.path.dirname(full_set_path)
    csv_out = os.path.join(out_dir, 'results_gpr.csv')
    map_out = os.path.join(out_dir, 'results_gpr.tif')
    fig_out = os.path.join(out_dir, 'results_gpr.png')
    
    assert os.path.exists(gpr_path), \
        f"Did not find `{gpr_path}` file.\nHave you trained the gaussian process regression (GPR) with lut.train_gpr()?"
    assert os.path.exists(full_set_path), \
        f"Did not find `{full_set_path}` file.\nHave you flattened the images with lut.image2csv()?"

    # --- 1. Load Model and Data ---
    # MATLAB: gpr = load(gpr_path); gprMdl = gpr.gprMdl;
    try:
        model_obj = joblib.load(gpr_path)
        gprMdl = model_obj['gprMdl']
        scaler = model_obj['scaler']
        feature_cols = model_obj['features']
        var_name = model_obj['var_name']
    except Exception as e:
        raise RuntimeError(f"Failed to load GPR model and/or scaler: {e}")

    val_in = pd.read_csv(full_set_path)
    
    # --- 2. Predict ---
    print("Working, usually takes around 1 minute")
    
    # Prepare input data: filter columns and scale
    X_val = val_in[feature_cols]
    X_val_scaled = scaler.transform(X_val)

    # MATLAB: res = predict(gprMdl, val_in);
    res, _ = gprMdl.predict(X_val_scaled, return_std=True)
    
    # --- 3. Save CSV Result ---
    # MATLAB: tab = array2table(res, 'VariableNames', {'Actot'});
    tab = pd.DataFrame(res, columns=[var_name])
    tab.to_csv(csv_out, index=False)
    print(f"saved `{csv_out}`")
    
    # --- 4. Reshape to Image and Save GeoTIFF ---
    # MATLAB: im = lut.csv2image_plane(val_in, res);
    im = csv2image_plane(val_in, res)
    
    # MATLAB: lut.write_tiff(im, map_out)
    write_tiff(im, map_out)
    
    # --- 5. Plot Image ---
    # MATLAB: lut.plot_image(im, 'GPR', fig_out)
    plot_image(im, 'GPR', fig_out)
    
if __name__ == '__main__':
    use_gpr()