import pandas as pd
import numpy as np
import os
from scipy.stats import qmc # MATLAB's lhsdesign equivalent

def latin_hypercube_input(tab=None, n_spectra=30, outdir=None):
    """
    Generates Latin Hypercube Sampling (LHS) for SCOPE input parameters.

    :param tab: DataFrame/Table with input borders (variable, lower, upper, include).
    :param n_spectra: Number of spectra/sets of parameters to generate.
    :param outdir: Output directory path.
    :return: None, saves to CSV.
    """
    
    if tab is None:
        outdir = os.path.join('input', 'dataset_for_verification')
        # Assuming the CSV path in MATLAB is a standard CSV/text file
        input_borders_path = os.path.join(outdir, 'input_borders.csv')
        tab = pd.read_csv(input_borders_path)
        
    if outdir is None:
        outdir = os.path.join('input', 'dataset_for_verification')
        
    out_file = os.path.join(outdir, 'lh_ts.csv')
    
    if os.path.exists(out_file):
        raise AssertionError(f"'{out_file}' file already exists, delete it first")
    
    # Filter by 'include' column (logical indexing)
    include = tab['include'].astype(bool)
    lb = tab[include]['lower'].values
    ub = tab[include]['upper'].values
    varnames = tab[include]['variable'].values
    
    n_vars = len(varnames)

    # Latin Hypercube Sampling
    # MATLAB's lhsdesign is equivalent to qmc.LatinHypercube in scipy.stats
    sampler = qmc.LatinHypercube(d=n_vars)
    lh = sampler.random(n=n_spectra)

    # Scale the samples to the bounds: params = (ub - lb) * lh + lb
    # The bounds need to be broadcastable with the lh matrix (n_spectra, n_vars)
    params = (ub - lb) * lh + lb

    if 'LIDFa' in varnames and 'LIDFb' in varnames:
        # abs(LIDFa + LIDFb) <= 1
        i_lidfa = np.where(varnames == 'LIDFa')[0][0]
        i_lidfb = np.where(varnames == 'LIDFb')[0][0]
        
        lidfa = params[:, i_lidfa]
        lidfb = params[:, i_lidfb]
        
        # MATLAB conversion logic:
        params[:, i_lidfa] = (lidfa + lidfb) / 2
        params[:, i_lidfb] = (lidfa - lidfb) / 2
        
    # Convert to DataFrame
    t = pd.DataFrame(params, columns=varnames)
    
    # Add time step 't' column
    t['t'] = np.arange(1, len(t) + 1)
    
    # Save to CSV
    t.to_csv(out_file, index=False)
    
    # Output message
    varnames_in = ', '.join(varnames)
    print(f"Sampled {len(varnames)} parameters: {varnames_in}")
    print(f"Saved lut input (parameters) in '{out_file}'")