import pandas as pd
import numpy as np
import os
import warnings
# For reading GeoTIFF/NetCDF (gdal/netcdf-based)
try:
    import rasterio
    import netCDF4
except ImportError:
    warnings.warn("Rasterio/NetCDF4 not installed. GeoTIFF/NetCDF reading will be mocked.")
    rasterio = None
    netCDF4 = None

# Mock functions for GeoTIFF/NetCDF reading if libraries are missing
def imread_single_band(path, var_name):
    """Mock/Wrapper for reading a single band."""
    if rasterio and path.endswith('.tif'):
        with rasterio.open(path) as src:
            return src.read(1)
    elif netCDF4 and path.endswith('.nc'):
        with netCDF4.Dataset(path, 'r') as nc_file:
            return nc_file.variables[var_name][:]
    else:
        # Mock behavior
        warnings.warn("Using mock image data.")
        return np.random.rand(100, 100) * 100

def image2csv(im_path=None, outdir=None):
    """
    // src/+lut/image2csv.m
    Reads image files (GeoTIFF or NetCDF), flattens them, cleans NaN pixels, 
    and saves the non-NaN pixels to a CSV file.
    """
    var_names = ['Cab', 'LAI']  # Default variables to process
    
    if im_path is None:
        im_path = os.path.join('..', 'exercise', 'images')
        outdir = os.path.join('..', 'exercise')
        name = 'full_set'
    else:
        name = os.path.splitext(os.path.basename(im_path))[0]
        
    csv_out = os.path.join(outdir, f'{name}.csv')
    
    # Determine if im_path is a directory or a single NetCDF file
    if os.path.isdir(im_path):
        assert os.path.exists(im_path), f"`{im_path}` does not exist. Please, put images (Cab.tif and LAI.tif) into that folder"
        ext = '.tif'
        get_val = lambda x: imread_single_band(os.path.join(im_path, f'{x}{ext}'), x)
    elif im_path.endswith('.nc'):
        ext = '.nc'
        get_val = lambda x: imread_single_band(im_path, x)
    else:
        raise ValueError(f"Unsupported image path format: {im_path}")

    print(f"reading images from `{im_path}`")
    
    # Get dimensions from the first variable
    try:
        v0 = get_val(var_names[0])
        r, c = v0.shape
    except Exception as e:
        raise RuntimeError(f"Could not read first image variable '{var_names[0]}': {e}")

    # 1. Flattening
    print(f"Flattening {r} rows, {c} columns")
    vals = np.full((r * c, len(var_names)), np.nan, dtype=np.float32)
    
    for i, var in enumerate(var_names):
        v = get_val(var)
        # MATLAB: vals(:, i) = v(:);  # flattening (column-major order is preserved here by Python's v.flat[:])
        vals[:, i] = v.flat[:]

    # 2. Clean NaN pixels
    i_nans = np.any(np.isnan(vals), axis=1)
    df = pd.DataFrame(vals, columns=var_names)
    df_clean = df[~i_nans].copy()
    
    # 3. Add pixel index column (1-based index)
    # MATLAB: df_clean.(['ind_' num2str(r), '_', num2str(c)]) = find(~i_nans);
    pixel_indices_1based = np.where(~i_nans)[0] + 1
    df_clean[f'ind_{r}_{c}'] = pixel_indices_1based
    
    print(f"found {sum(~i_nans)} not-nan pixels")
    
    # 4. Saving
    df_clean.to_csv(csv_out, index=False)
    print(f"saved to `{csv_out}`")

if __name__ == '__main__':
    image2csv()