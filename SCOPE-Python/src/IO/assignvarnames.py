def assignvarnames():
    """
    Translation of src/IO/assignvarnames.m
    Initializes a list of dictionaries (equivalent to a MATLAB struct array)
    to store variable names for SCOPE parameters.
    """
    # Initialize V, analogous to: V = struct('Name','','Val', zeros(64,1))
    V = []

    # Leaf Biochemical/Optical Parameters (1-11)
    V.append({'Name': 'Cab', 'Val': 0.0})
    V.append({'Name': 'Cca', 'Val': 0.0})
    V.append({'Name': 'Cdm', 'Val': 0.0})
    V.append({'Name': 'Cw', 'Val': 0.0})
    V.append({'Name': 'Cs', 'Val': 0.0})
    V.append({'Name': 'N', 'Val': 0.0})
    V.append({'Name': 'Cant', 'Val': 0.0})       # Added March 2017
    V.append({'Name': 'Cp', 'Val': 0.0})
    V.append({'Name': 'Cbc', 'Val': 0.0})       # carbon based consituents PROSPECT-PRO
    V.append({'Name': 'rho_thermal', 'Val': 0.0})
    V.append({'Name': 'tau_thermal', 'Val': 0.0})

    # Biochemical/Physiological Parameters (12-21)
    V.append({'Name': 'Vcmax25', 'Val': 0.0})
    V.append({'Name': 'BallBerrySlope', 'Val': 0.0})
    V.append({'Name': 'BallBerry0', 'Val': 0.0})
    V.append({'Name': 'Type', 'Val': ''})       # Note: This is stored as a string ('C3'/'C4') in the code
    V.append({'Name': 'kV', 'Val': 0.0})
    V.append({'Name': 'Rdparam', 'Val': 0.0})
    V.append({'Name': 'fqe', 'Val': 0.0})
    V.append({'Name': 'Kn0', 'Val': 0.0})
    V.append({'Name': 'Knalpha', 'Val': 0.0})
    V.append({'Name': 'Knbeta', 'Val': 0.0})

    # Canopy Parameters (22-28)
    V.append({'Name': 'LAI', 'Val': 0.0})
    V.append({'Name': 'hc', 'Val': 0.0})
    V.append({'Name': 'zo', 'Val': 0.0})
    V.append({'Name': 'd', 'Val': 0.0})
    V.append({'Name': 'LIDFa', 'Val': 0.0})
    V.append({'Name': 'LIDFb', 'Val': 0.0})
    V.append({'Name': 'leafwidth', 'Val': 0.0})

    # Meteorological Inputs (29-37)
    V.append({'Name': 'z', 'Val': 0.0})
    V.append({'Name': 'Rin', 'Val': 0.0})
    V.append({'Name': 'Ta', 'Val': 0.0})
    V.append({'Name': 'Rli', 'Val': 0.0})
    V.append({'Name': 'p', 'Val': 0.0})
    V.append({'Name': 'ea', 'Val': 0.0})
    V.append({'Name': 'u', 'Val': 0.0})
    V.append({'Name': 'Ca', 'Val': 0.0})
    V.append({'Name': 'Oa', 'Val': 0.0})

    # Canopy/Soil Resistances and Constants (38-45)
    V.append({'Name': 'rb', 'Val': 0.0})
    V.append({'Name': 'Cd', 'Val': 0.0})
    V.append({'Name': 'CR', 'Val': 0.0})
    V.append({'Name': 'CD1', 'Val': 0.0})
    V.append({'Name': 'Psicor', 'Val': 0.0})
    V.append({'Name': 'CSSOIL', 'Val': 0.0})
    V.append({'Name': 'rbs', 'Val': 0.0})
    V.append({'Name': 'rwc', 'Val': 0.0})

    # Time/Site/Angle/Soil Parameters (46-53)
    V.append({'Name': 'startDOY', 'Val': 0.0})
    V.append({'Name': 'endDOY', 'Val': 0.0})
    V.append({'Name': 'LAT', 'Val': 0.0})
    V.append({'Name': 'LON', 'Val': 0.0})
    V.append({'Name': 'timezn', 'Val': 0.0})
    V.append({'Name': 'tts', 'Val': 0.0})
    V.append({'Name': 'tto', 'Val': 0.0})
    V.append({'Name': 'psi', 'Val': 0.0})

    # Soil/Stress/Model Parameters (54-59)
    V.append({'Name': 'SMC', 'Val': 0.0})
    V.append({'Name': 'Tyear', 'Val': 0.0})
    V.append({'Name': 'beta', 'Val': 0.0})
    V.append({'Name': 'kNPQs', 'Val': 0.0})
    V.append({'Name': 'qLs', 'Val': 0.0})
    V.append({'Name': 'stressfactor', 'Val': 0.0})

    # BSM/Spectral (60-63)
    V.append({'Name': 'spectrum', 'Val': 0.0})
    V.append({'Name': 'BSMBrightness', 'Val': 0.0})
    V.append({'Name': 'BSMlat', 'Val': 0.0})
    V.append({'Name': 'BSMlon', 'Val': 0.0})

    # Soil Thermal/Moisture (64-68)
    V.append({'Name': 'rss', 'Val': 0.0})
    V.append({'Name': 'rs_thermal', 'Val': 0.0})
    V.append({'Name': 'cs', 'Val': 0.0})
    V.append({'Name': 'rhos', 'Val': 0.0})
    V.append({'Name': 'lambdas', 'Val': 0.0})

    # The MATLAB file initializes 64 entries but lists 68 named entries.
    # The Python list accommodates all 68 named entries.
    return V