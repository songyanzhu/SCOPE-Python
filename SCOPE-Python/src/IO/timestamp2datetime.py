import numpy as np
import pandas as pd
import warnings
import re
from datetime import datetime

def timestamp2datetime(ts, year_n=None):
    """
    Translation of src/IO/timestamp2datetime.m
    Converts various timestamp formats (DOY, yyyyMMddHHmmss...) to datetime objects.
    """
    
    # 1. Convert to NumPy array and handle NaNs
    if not isinstance(ts, np.ndarray):
        ts = np.array(ts, dtype=object)
    
    # Handle MATLAB's -9999 as NaN/NaT
    if ts.dtype == 'object':
        ts = np.where(ts == '-9999', np.nan, ts)
        ts = np.where(ts == -9999, np.nan, ts)
        ts = np.where(pd.isna(ts), pd.NaT, ts)
    elif np.issubdtype(ts.dtype, np.number):
        ts = np.where(ts == -9999, np.nan, ts)
        ts = np.where(np.isnan(ts), pd.NaT, ts)

    # Convert to Pandas Series of strings for uniform format checking
    ts_str = pd.Series(ts).astype(str).str.replace('<NA>', '').replace('nan', '')
    
    # 2. DOY Conversion (if all values look like DOY)
    if year_n is not None and not pd.isna(ts_str).all() and np.all([re.fullmatch(r'\d+', s) for s in ts_str if s]) and np.max([int(s) for s in ts_str.dropna().unique() if s]) <= 367:
        warnings.warn(f't is DOY, converting to date with year = {year_n[0] if year_n.ndim > 0 else year_n}')
        
        dt_out = []
        for d, y in zip(ts, year_n):
            if pd.isna(d) or pd.isna(y):
                dt_out.append(pd.NaT)
            else:
                try:
                    # pd.to_datetime(f'{year}').toordinal() + doy
                    dt_out.append(pd.to_datetime(f'{int(y)}', format='%Y') + pd.to_timedelta(d - 1, unit='D'))
                except ValueError:
                    dt_out.append(pd.NaT)
        return np.array(dt_out)

    # 3. Determine Format and Convert
    # Find the length of the longest non-empty timestamp string
    valid_ts = ts_str[ts_str.str.len() > 0]
    if valid_ts.empty:
        return np.full_like(ts, pd.NaT)

    ts_len = valid_ts.str.len().max()
    
    # MATLAB: switch size(ts, 2)
    format_map = {
        18: '%Y%m%d%H%M%S.%f', # SS.SSS
        14: '%Y%m%d%H%M%S',  # SS
        12: '%Y%m%d%H%M',    # MM
        10: '%Y%m%d%H',      # HH (non-standard)
        8:  '%Y%m%d',      # DD
    }
    
    # Find the closest matching format length
    format_key = min(format_map.keys(), key=lambda x: abs(x - ts_len))
    if format_key not in format_map:
        raise ValueError('Format of timestamp is not any truncation of `yyyyMMddHHmmss`')

    # Convert to datetime using the detected format
    try:
        dt = pd.to_datetime(ts_str, format=format_map[format_key], errors='coerce').values
    except ValueError as e:
        warnings.warn(f"Failed to parse timestamp with format {format_map[format_key]}: {e}")
        dt = np.full_like(ts, pd.NaT) # Return NaT array

    # Check for NaT resulting from bad conversion
    if np.any(pd.isna(dt)) and format_key == 18:
        warnings.warn('milliseconds timestamp is not perfectly matched and resulted in NaTs.')
        
    return dt