import numpy as np
import pandas as pd

# --- Mock/Helper Functions for I/O ---

def calczenithangle(DOY, time_hr, LAT, LON, LON_site, LAT_site):
    """
    MOCK: Placeholder for solar zenith angle calculation.
    Returns a constant 45 degrees (pi/4 radians).
    """
    return np.full_like(DOY, np.pi/4)

def satvap(Ta_celsius):
    """
    Calculates saturated vapor pressure using Magnus-Tetens formula (hPa).
    """
    return 6.107 * np.power(10, 7.5 * Ta_celsius / (237.3 + Ta_celsius))

def juliandate(dt):
    """
    Converts datetime/timestamp to Julian date number (ordinal + 1).
    Assumes dt is a Pandas/NumPy datetime array.
    """
    if isinstance(dt, pd.Timestamp):
        return dt.toordinal() + 1
    elif isinstance(dt, np.datetime64):
        # Convert to Pandas for easy ordinal
        return pd.to_datetime(dt).toordinal() + 1
    else:
        # Assume it's an array of dates
        return pd.to_datetime(dt).toordinal().values + 1