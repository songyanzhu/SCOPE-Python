import pandas as pd
import numpy as np
import re
import warnings

def csv2image_plane(val_in, res):
    """
    // src/+lut/csv2image_plane.m
    Reshapes flattened CSV data back into an image plane (2D array).
    """
    
    # MATLAB: val_in is a table/DataFrame
    vars_names = val_in.columns
    
    # 1. Find the index column name (e.g., 'ind_R_C')
    # MATLAB: i_ind = ~cellfun(@isempty, strfind(vars, 'ind')); 
    ind_col_name = None
    for name in vars_names:
        if 'ind' in name:
            ind_col_name = name
            break
            
    assert ind_col_name is not None, "Input table must contain an index column (e.g., 'ind_R_C')."
    
    # 2. Extract dimensions R and C from the column name
    # MATLAB: splt = strsplit(vars{i_ind}, '_'); r = str2num(splt{2}); c = str2num(splt{3});
    match = re.search(r'ind_(\d+)_(\d+)', ind_col_name)
    if not match:
        # Fallback if the naming convention is slightly off, but 'ind' exists
        warnings.warn("Could not parse image dimensions (R, C) from index column name. Using largest pixel index to estimate R*C.")
        pix_ind_max = val_in[ind_col_name].max()
        r, c = int(np.sqrt(pix_ind_max)), int(pix_ind_max / np.sqrt(pix_ind_max)) # Crude guess
    else:
        r = int(match.group(1))
        c = int(match.group(2))
        
    # 3. Create NaN image array
    im = np.full((r, c), np.nan, dtype=np.float32)
    
    # 4. Fill image array at non-NaN pixel indices
    # MATLAB: pix_ind = table2array(val_in(:, i_ind));
    # MATLAB uses 1-based indexing for im(pix_ind). Python uses 0-based.
    pix_ind_1based = val_in[ind_col_name].values
    pix_ind_0based = pix_ind_1based - 1 
    
    # Check if 'res' is a DataFrame/Series or a plain array. It should be an array of values.
    if isinstance(res, (pd.Series, pd.DataFrame)):
        # Assuming res is a single column/series
        res = res.squeeze().values

    # Res must be a vector of the same length as the number of non-NaN pixels
    assert len(res) == len(pix_ind_0based), "Result vector 'res' length does not match pixel index count."
    
    im.flat[pix_ind_0based] = res
    
    return im