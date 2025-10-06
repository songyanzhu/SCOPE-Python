import numpy as np
import pandas as pd
import os
from scipy.interpolate import interp1d
import warnings

# Assume helper functions are available
from .timestamp2datetime import timestamp2datetime
from .load_mSCOPE_ts import load_mSCOPE_ts
from .timestamp_utils import calczenithangle, satvap # Mock/Placeholder for these

# --- Placeholder/Mock Functions ---

# Placeholder for solar angle calculation (calczenithangle)
def calczenithangle(DOY, time_hr, latitude, longitude, lon, lat):
    """MOCK: Returns a fixed zenith angle in radians (e.g., 45 deg = pi/4)"""
    return np.full_like(DOY, np.pi/4)

# Placeholder for saturation vapor pressure (satvap)
def satvap(Ta_celsius):
    """MOCK: Returns saturation vapor pressure in hPa"""
    return 6.107 * np.power(10, 7.5 * Ta_celsius / (237.3 + Ta_celsius))

# --- Main Function ---

def load_timeseries(V, F, xyt, path_input):
    """
    Translation of src/IO/load_timeseries.m
    Loads, filters, and prepares all time-series input data.
    """
    
    # --- Filenames and Setup ---
    F_dict = {f['FileID']: f['FileName'] for f in F if 'FileID' in f and 'FileName' in f}
    
    Dataset_dir = f"dataset {F_dict.get('Dataset_Dir', '')}"
    meteo_ec_csv = F_dict.get('meteo_ec_csv')
    vegetation_retrieved_csv = F_dict.get('vegetation_retrieved_csv')
    mSCOPE_csv = F_dict.get('mSCOPE_csv')

    t_column = F_dict.get('t')
    year_column = F_dict.get('year')
    atmo_column = F_dict.get('atmos_names')
    
    # --- Read main dataset ---
    if not meteo_ec_csv:
        raise ValueError("meteo_ec_csv filename must be provided in input structure F.")

    try:
        df = pd.read_csv(os.path.join(path_input, Dataset_dir, meteo_ec_csv), 
                          na_values=['.', 'NA', 'N/A', '-9999'])
    except Exception as e:
        raise IOError(f"Failed to read meteo/EC file: {e}")

    # --- Time parsing and standardization ---
    
    # MATLAB: t_ = df.(t_column)
    t_ = df[t_column].values
    
    # Convert DOY to datetime if necessary
    if np.all(t_ <= 367) and t_.ndim > 0 and t_.size > 0:
        if not year_column or df[year_column].isnull().all():
            warnings.warn('t is DOY, converting to date with year = 2020, as `year` in .csv was empty')
            year_n_col = np.full_like(t_, 2020)
        else:
            year_n_col = df[year_column].values
        
        # MATLAB: t_ = datestr(datenum(year_n, 0, t_), 'yyyymmddHHMMSS.FFF');
        # t_ now becomes a list of datetime objects (or comparable numerical representation)
        t_dt = []
        for doy, year in zip(t_, year_n_col):
            if np.isnan(doy) or np.isnan(year):
                t_dt.append(pd.NaT)
            else:
                try:
                    t_dt.append(pd.to_datetime(f'{int(year)}', format='%Y') + pd.to_timedelta(doy - 1, unit='D'))
                except ValueError:
                    t_dt.append(pd.NaT)
        t_ = np.array(t_dt)
    else:
        t_ = timestamp2datetime(t_)

    # Convert start/end DOY to datetime for filtering
    xyt['startDOY'] = timestamp2datetime(xyt['startDOY'])
    xyt['endDOY'] = timestamp2datetime(xyt['endDOY'])
    
    # --- Filtering ---
    if np.any(pd.isna(t_)):
         raise ValueError('Some time series values contain NaN/NaT. Please remove them as the date must be known.')
         
    # MATLAB: time_i = (t_ >= xyt.startDOY) & (t_ <= xyt.endDOY);
    # NumPy/Pandas handles the datetime comparison
    time_i = (t_ >= xyt['startDOY']) & (t_ <= xyt['endDOY'])
    
    df_sub = df[time_i].reset_index(drop=True)
    t_ = t_[time_i]
    year_n = pd.to_datetime(t_).year.values # Ensure we keep the year array for legacy

    xyt['t'] = t_
    xyt['year'] = year_n

    # --- Optional Interpolation/mSCOPE files ---
    
    mly_ts = {}
    interpolatable_cols = {}
    
    if vegetation_retrieved_csv:
        df_int = pd.read_csv(os.path.join(path_input, Dataset_dir, vegetation_retrieved_csv), 
                              na_values=['.', 'NA', 'N/A', '-9999'])
        t_int = timestamp2datetime(df_int[t_column].values)
        
        # Check time bounds
        if t_int.min() > t_.max() or t_int.max() < t_.min():
            raise ValueError('Correct the timestamp. vegetation_retrieved_csv time range is outside meteo_ec_csv range.')
        if t_int.min() > t_.min():
            warnings.warn('vegetation_retrieved_csv starts later, head simulations will be 0 or NaNs')
        if t_int.max() < t_.max():
            warnings.warn('vegetation_retrieved_csv ends earlier, tail simulations will be 0 or NaNs')

        t_int_sec = (t_int - t_int.min()).total_seconds().values
        t_sec = (t_ - t_int.min()).total_seconds().values
        
        # Prepare interpolation functions (store function in dict for later use)
        for col in df_int.columns.drop(t_column):
            valid_idx = df_int[col].notna() & np.isfinite(t_int_sec)
            if valid_idx.sum() > 1:
                interpolatable_cols[col] = interp1d(t_int_sec[valid_idx], df_int[col].values[valid_idx], 
                                                    kind='linear', fill_value=np.nan, bounds_error=False)
            else:
                 interpolatable_cols[col] = lambda x: np.nan # No valid data for interpolation

    if mSCOPE_csv:
        mSCOPE_ts_path = os.path.join(path_input, Dataset_dir, mSCOPE_csv)
        nly = F_dict.get('nly', 1) # Assumes F(11).FileName contains nly
        if isinstance(nly, str): nly = int(nly)
        
        mly_ts = load_mSCOPE_ts(mSCOPE_ts_path, nly, t_column, t_)

    # --- Read fields and assign to V ---
    
    v_names = [v['Name'] for v in V]
    atmo_paths = []
    
    for i, v_name in enumerate(v_names):
        if v_name in F_dict:
            col_name = F_dict[v_name]
            
            if col_name == atmo_column:
                # Handle Irradiance files separately
                atmo_paths = os.path.join(path_input, Dataset_dir, df_sub[col_name].values)
                continue
            
            # Interpolation from vegetation file
            if col_name in interpolatable_cols:
                V[i]['Val'] = interpolatable_cols[col_name](t_sec)
            # Direct read from meteo/EC file
            else:
                tmp = df_sub[col_name].values
                if np.all(np.isnan(tmp)):
                    warnings.warn(f"{col_name} has NaNs along all timestamp. Calculations may fail")
                V[i]['Val'] = tmp

    # --- Special Cases and Conversions ---
    
    # 1. tts calculation (Sun Zenith Angle)
    if 'tts' not in F_dict:
        vi_tts = v_names.index('tts')
        DOY_ = (t_ - pd.to_datetime(t_.astype(str)).floor('D')).dt.days.values + 1 # Day of Year (1-based)
        time_ = (t_ - pd.to_datetime(t_.astype(str)).floor('D')).dt.total_seconds().values / 3600 # Hours since midnight
        
        # Fallback for midday time if all times are 0
        if np.all(time_ == 0):
            time_ = np.full_like(time_, 12.0)
            warnings.warn('Taking midday time to calculate tts')
        
        lon, lat, timezn = xyt['LON'], xyt['LAT'], xyt['timezn']
        
        # ttsR: sun zenith angle in rad
        ttsR = calczenithangle(DOY_, time_ - timezn, lat, lon, lon, lat) # Arguments match MATLAB's expected
        
        # Convert to degrees, max 85
        V[vi_tts]['Val'] = np.minimum(85, ttsR / np.pi * 180)

    # 2. ea calculation (Vapor Pressure)
    if 'ea' not in F_dict and 'Ta' in F_dict:
        vi_ea = v_names.index('ea')
        vi_ta = v_names.index('Ta')
        ta = V[vi_ta]['Val']
        es = satvap(ta)
        
        rh_column = F_dict.get('RH')
        vpd_column = F_dict.get('VPD')
        
        if rh_column:
            rh = df_sub[rh_column].values
            rh = np.where(rh == -9999, np.nan, rh)
            if np.any(rh > 10):
                rh /= 100
                warnings.warn('Converted relative humidity from [0 100] to [0 1]')
            V[vi_ea]['Val'] = es * rh
            warnings.warn('Calculated ea from Ta and RH')
        elif vpd_column:
            vpd = df_sub[vpd_column].values
            vpd = np.where(vpd == -9999, np.nan, vpd)
            V[vi_ea]['Val'] = es - vpd
            warnings.warn('Calculated ea from Ta and VPD')
            if np.any(V[vi_ea]['Val'] < 0):
                warnings.warn('Some ea < 0, is your VPD in hPa?')

    # 3. Unit Conversions
    # p (Pressure)
    if 'p' in F_dict:
        vi_p = v_names.index('p')
        if np.any(V[vi_p]['Val'] < 500):
            V[vi_p]['Val'] *= 10
            warnings.warn('Converted air pressure from kPa to hPa')

    # SMC (Soil Moisture Content)
    if 'SMC' in F_dict:
        vi_smc = v_names.index('SMC')
        if np.any(V[vi_smc]['Val'] > 1):
            V[vi_smc]['Val'] /= 100
            warnings.warn('Converted soil moisture content from from [0 100] to [0 1]')

    return V, xyt, mly_ts, atmo_paths