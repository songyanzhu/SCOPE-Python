import pandas as pd
import numpy as np
import time
import warnings

def _top_indices(i, meas_arr, lut_arr):
    """
    // src/+lut/lut_search.m -> inner function top_indices(i, meas_arr, lut_arr)
    Finds the indices of the top 10 best matching LUT entries (smallest RMSE).
    """
    # Note: 'i' is the 1-based index of the measurement pixel being queried
    # We use 0-based indexing for meas_arr
    meas_row = meas_arr[i - 1, :]
    
    # Calculate RMSE (Euclidean distance)
    # MATLAB: rmses = sqrt(mean((meas_arr(i, :) - lut_arr) .^ 2, 2));
    rmses = np.sqrt(np.mean((meas_row - lut_arr) ** 2, axis=1))
    
    # Sort and get top 10 indices (1-based to match MATLAB's use with response array)
    # MATLAB: [~, I] = sort(rmses); ind = I(1:10)';
    I = np.argsort(rmses)
    
    # Select top 10 indices (I is 0-based, so add 1 for MATLAB-style indexing into response)
    ind_0based = I[:10]
    
    # Convert to 1-based for consistent indexing if 'response' array is 1-based
    # If 'response' is 0-based (Python standard), we must return 0-based indices.
    # We return 0-based indices and assume 'response' array is 0-based.
    return ind_0based

def lut_search(params, lut_params, response):
    """
    // src/+lut/lut_search.m
    Performs the LUT search using the smallest Root Mean Square Error (RMSE) 
    between input parameters and LUT parameters.
    """
    
    p_names = params.columns.tolist()
    lut_names = lut_params.columns.tolist()
    
    # --- 1. Normalization and Alignment ---
    
    # Create the normalized input array aligned with LUT parameter names
    p_ordered = np.full((params.shape[0], len(lut_names)), np.nan, dtype=np.float32)

    lut_params_norm = lut_params.copy()
    
    # Note: This loop performs min/max calculation and normalization on both tables
    for i, v in enumerate(lut_names):
        if v not in params.columns:
            warnings.warn(f"Parameter '{v}' is in LUT but not in input. Skipping for search.", UserWarning)
            continue
            
        v_min = lut_params[v].min()
        v_max = lut_params[v].max()
        
        # Normalize LUT
        lut_params_norm[v] = (lut_params_norm[v] - v_min) / (v_max - v_min)
        
        # Normalize Input (p_ordered)
        p_ordered[:, i] = (params[v].values - v_min) / (v_max - v_min)
        
    lut_arr = lut_params_norm.values
    
    # --- 2. LUT Search (Vectorized Query) ---

    print("Starting LUT search...")
    start_time = time.time()
    
    n_pixels = p_ordered.shape[0]
    
    # Query each pixel in the input (params) against the entire LUT (lut_arr)
    # MATLAB: ind = arrayfun(@(x) top_indices(x, p_ordered, lut_arr), 1:size(p_ordered, 1), 'UniformOutput', false);
    
    # Query: 1-based index (i) to 0-based index (i-1)
    indices = []
    for i in range(1, n_pixels + 1):
        if i % 100000 == 0:
            print(f"done {i} / {n_pixels} pixels")
            
        # _top_indices returns 0-based indices into the LUT (response array)
        ind_i = _top_indices(i, p_ordered, lut_arr)
        indices.append(ind_i)
        
    ind = np.array(indices)
    
    end_time = time.time()
    print(f"LUT search finished in {end_time - start_time:.2f} seconds.")

    # --- 3. Result Aggregation ---
    
    # response is the vector of SCOPE outputs from the LUT (e.g., Actot)
    # MATLAB: actot = response(ind);
    # Since 'ind' is 0-based, we can use it to directly index the 'response' array
    actot = response[ind]
    
    # MATLAB: res = nanmedian(actot, 2); (Median across the top 10 matches)
    res = np.nanmedian(actot, axis=1)
    
    # MATLAB: res_std = nanstd(actot, 0, 2);
    res_std = np.nanstd(actot, axis=1)

    return res, res_std