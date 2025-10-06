import pandas as pd
import numpy as np
import os
import warnings

# Mock external functions
from .lut_search import lut_search
from .csv2image_plane import csv2image_plane
from .plot_image import plot_image

# Assumed external function to write image
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

def use_rmse(lut_in_path=None, lut_out_path=None, full_set_path=None, var_name='Actot'):
    """
    // src/+lut/use_rmse.m
    Performs LUT search (RMSE matching) on a full dataset.
    """
    if lut_in_path is None:
        in_dir = os.path.join('..', 'exercise')
        lut_in_path = os.path.join(in_dir, 'lut_in.csv')
        lut_out_path = os.path.join(in_dir, 'lut_out.csv')
        full_set_path = os.path.join(in_dir, 'full_set.csv')
        
    out_dir = os.path.dirname(full_set_path)
    csv_out = os.path.join(out_dir, 'results_rmse.csv')
    map_out = os.path.join(out_dir, 'results_rmse.tif')
    fig_out = os.path.join(out_dir, 'results_rmse.png')
    
    assert os.path.exists(lut_in_path), \
        f"Did not find `{lut_in_path}` file.\nHave you generated the LUT input with lut.generate_lut_input()?"
    assert os.path.exists(lut_out_path), \
        f"Did not find `{lut_out_path}` file.\nHave you copied the `fluxes.csv` from `../output` after SCOPE run on lut_in.csv?"
    assert os.path.exists(full_set_path), \
        f"Did not find `{full_set_path}` file.\nHave you flattened the images with lut.image2csv()?"

    # --- 1. Load Data ---
    lut_in = pd.read_csv(lut_in_path)
    
    # MATLAB: opt = detectImportOptions(lut_out_path); flu = readtable(lut_out_path, opt);
    flu = pd.read_csv(lut_out_path)
    
    # MATLAB: lut_out = flu.(var_name);
    if var_name not in flu.columns:
        raise ValueError(f"Variable '{var_name}' not found in LUT output file.")
    lut_out = flu[var_name].values
    
    val_in = pd.read_csv(full_set_path)
    
    # --- 2. LUT Search ---
    # [res, res_std] = lut.lut_search(val_in, lut_in, lut_out);
    res, res_std = lut_search(val_in, lut_in, lut_out)
    
    # --- 3. Save CSV Result ---
    # MATLAB: tab = array2table([res, res_std], 'VariableNames', {'Actot', 'Actot_sd'});
    tab = pd.DataFrame({'Actot': res, 'Actot_sd': res_std})
    tab.to_csv(csv_out, index=False)
    print(f"saved `{csv_out}`")
    
    # --- 4. Reshape to Image and Save GeoTIFF ---
    # MATLAB: im = lut.csv2image_plane(val_in, res);
    im = csv2image_plane(val_in, res)
    
    # MATLAB: lut.write_tiff(im, map_out)
    write_tiff(im, map_out)
    
    # --- 5. Plot Image ---
    # MATLAB: lut.plot_image(im, 'RMSE', fig_out)
    plot_image(im, 'RMSE', fig_out)

if __name__ == '__main__':
    use_rmse()