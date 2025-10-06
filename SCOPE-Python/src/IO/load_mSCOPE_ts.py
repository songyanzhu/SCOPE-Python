import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import warnings

# Assume timestamp2datetime is available globally or imported
from .timestamp2datetime import timestamp2datetime

def load_mSCOPE_ts(mSCOPE_ts_path, nly, t_column, t_):
    """
    Translation of src/IO/load_mSCOPE_ts.m
    Loads and interpolates multilayer time series data.
    """
    
    try:
        # MATLAB: mly_df = readtable(mSCOPE_ts_path);
        # We need to ensure we read missing values correctly, assuming NaN for missing.
        mly_df = pd.read_csv(mSCOPE_ts_path, na_values=['.', 'NA', 'N/A', '-9999'])
    except Exception as e:
        raise IOError(f"Failed to read mSCOPE time series file {mSCOPE_ts_path}: {e}")

    # --- Time conversion and validation ---
    
    t_mly = mly_df[t_column]
    t_ = np.array(t_) # Ensure target time is a NumPy array for comparison
    
    if t_mly.max() <= 367: # doy is provided
        # Requires 't_' to contain datetime objects to extract the year
        if not np.all([isinstance(t_obj, pd.Timestamp) for t_obj in t_]):
             warnings.warn("Target time 't_' does not appear to be datetimes; assuming year=2020 for DOY conversion.")
             year_n = 2020
        else:
             year_n = pd.to_datetime(t_).year.unique()
        
             if len(year_n) != 1:
                 raise ValueError('Multiple seasons in mSCOPE are supported only if t_column is TIMESTAMP, not DOY')
             year_n = year_n[0]
        
        # Convert DOY to datetime strings (MATLAB: datestr(datenum(year_n, 0, t_mly), 'yyyymmddHHMMSS.FFF'))
        t_mly_dt = [pd.to_datetime(f'{year_n}').toordinal() + doy for doy in t_mly.values]
        t_mly = pd.to_datetime(t_mly_dt, origin='toordinal', unit='D')
    else:
        # Assume t_column is already a TIMESTAMP string/number
        t_mly = timestamp2datetime(t_mly.values) # Ensure consistent datetime objects

    # Ensure all times are comparable (e.g., converted to seconds since epoch)
    t_mly_sec = (t_mly - t_mly.min()).total_seconds().values
    t_sec = (pd.to_datetime(t_) - t_mly.min()).total_seconds().values
    
    if t_sec.min() < t_mly_sec.min() or t_sec.max() > t_mly_sec.max():
        raise ValueError('Interpolation of mSCOPE values is not possible. Time range from mSCOPE file does not fully cover range of TS data.')
    
    # --- Interpolation ---
    
    # MATLAB: mly_df.(t_column) = []; (remove time column)
    mly_df = mly_df.drop(columns=[t_column])
    
    # MATLAB: variables = mly_df.Properties.VariableNames;
    variables = mly_df.columns.tolist()
    
    # MATLAB: param_names = variables(1:nly:length(variables));
    param_names = variables[::nly]
    
    # MATLAB: mly_in_t_ = interp1(t_mly, table2array(mly_df), t_);
    mly_in_t_ = np.zeros((len(t_sec), mly_df.shape[1]))
    
    # Interpolate each parameter column independently
    for i, col_name in enumerate(mly_df.columns):
        # interp1d requires a 1D input for the function values
        f_interp = interp1d(t_mly_sec, mly_df[col_name].values, kind='linear', fill_value=np.nan, bounds_error=False)
        mly_in_t_[:, i] = f_interp(t_sec)
        
    # --- Reshape and Structure Output ---

    # MATLAB: mly_split = mat2cell(mly_in_t_, length(t_), repmat(nly, 1, length(param_names)));
    # MATLAB: mly_ts = cell2struct(mly_split, param_names, 2);
    
    mly_ts = {}
    
    # The columns are [P1_L1, P1_L2, ..., P1_Ln, P2_L1, ...]
    for i, p_name in enumerate(param_names):
        start_col = i * nly
        end_col = (i + 1) * nly
        # mly_ts['p_name'] = [n_time_steps, nly] array
        mly_ts[p_name] = mly_in_t_[:, start_col:end_col]

    mly_ts['nly'] = nly
    
    return mly_ts