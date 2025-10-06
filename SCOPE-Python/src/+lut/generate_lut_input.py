import pandas as pd
import numpy as np
import os
import warnings
from datetime import datetime, timedelta
# For lhsdesign, use scipy.stats.qmc.LHS
try:
    from scipy.stats.qmc import LatinHypercube
except ImportError:
    warnings.warn("SciPy's QMC module (needed for lhsdesign) is not available.")
    LatinHypercube = None

def generate_lut_input(tab=None, n_spectra=None, outdir=None):
    """
    // src/+lut/generate_lut_input.m
    Generates input parameters for the LUT using Latin Hypercube Sampling (LHS).
    """
    if tab is None:
        # MATLAB: tab = readtable('+lut/input_borders.csv');
        # Assuming input_borders.csv is in the same directory as this module's parent ('lut')
        csv_path = os.path.join(os.path.dirname(__file__), 'input_borders.csv')
        # Use pandas to read the CSV
        tab = pd.read_csv(csv_path, sep='\t')
        n_spectra = 1000
        outdir = os.path.join('..', 'exercise')

    out_file = os.path.join(outdir, 'lut_in.csv')
    
    assert not os.path.exists(out_file), f"`{out_file}` file already exists, delete it first"
    
    # 1. Filter parameters to include
    include = tab['include'].astype(bool)
    lb = tab[include]['lower'].values
    ub = tab[include]['upper'].values
    varnames = tab[include]['variable'].tolist()

    if LatinHypercube is None:
        raise ImportError("LHS design requires scipy.stats.qmc. Please install SciPy.")

    # 2. Latin Hypercube Sampling
    # MATLAB: lh = lhsdesign(n_spectra, sum(include));
    # LHS generates samples scaled from 0 to 1
    # Note: lhsdesign is 0-centered random, LatinHypercube is quasi-random/deterministic
    lhs = LatinHypercube(d=len(lb), seed=42)
    lh = lhs.random(n=n_spectra)

    # 3. Rescale from [0, 1] to [lb, ub]
    # MATLAB: params = (ub-lb) .* lh + lb;
    params = (ub - lb) * lh + lb
    
    # 4. Handle LIDF parameters (LIDFa and LIDFb) constraint
    if 'LIDFa' in varnames:
        # abs(LIDFa + LIDFb) <= 1
        i_lidfa = varnames.index('LIDFa')
        i_lidfb = varnames.index('LIDFb')
        
        # MATLAB swaps the parameter roles to satisfy the constraint
        # new_LIDFa = (old_LIDFa + old_LIDFb) / 2
        # new_LIDFb = (old_LIDFa - old_LIDFb) / 2
        lidfa_old = params[:, i_lidfa]
        lidfb_old = params[:, i_lidfb]
        
        params[:, i_lidfa] = (lidfa_old + lidfb_old) / 2
        params[:, i_lidfb] = (lidfa_old - lidfb_old) / 2

    # 5. Create DataFrame and add time variable
    t = pd.DataFrame(params, columns=varnames)
    
    # Time variable generation (t)
    # MATLAB: t = datestr(datetime(2022, 7, 1, 0, 12, 1:n_spectra), 'yyyymmddHHMMSS.FFF');
    start_time = datetime(2022, 7, 1, 0, 12, 0)
    time_list = [start_time + timedelta(seconds=i) for i in range(1, n_spectra + 1)]
    # Python datetime formatting: 'yyyymmddHHMMSS.FFF' -> '%Y%m%d%H%M%S.%f'
    t['t'] = [dt.strftime('%Y%m%d%H%M%S.%f')[:17] for dt in time_list]

    # 6. Save to CSV
    # MATLAB: writetable(t, out_file)
    t.to_csv(out_file, index=False)
    
    varnames_in = ', '.join(varnames)
    print(f"Sampled {len(varnames)} parameters: {varnames_in}")
    print(f"Saved lut input (parameters) in `{out_file}`")

if __name__ == '__main__':
    generate_lut_input()