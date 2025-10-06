import numpy as np

def initialize_output_structures(spectral):
    """
    Translation of src/IO/initialize_output_structures.m
    Initializes the main output structures (rad, thermal, fluxes) with NaN.
    """

    # Initialization of scalar fields to NaN
    scalar_fields = [
        'Rntot', 'lEtot', 'Htot', 'Atot', 'Rnctot', 'lEctot', 'Hctot', 'Actot',
        'Rnstot', 'lEstot', 'Hstot', 'Gtot', 'Resp', 'Tcave', 'Tsave', 
        'raa', 'rawc', 'raws', 'ustar', 'Lout', 'Loutt', 'Eoutte', 'PAR', 
        'Eouto', 'Eout', 'Lout', 'Lo_', 'Ta' # MATLAB also sets Lo_ and Lo_ (incorrectly) to single NaN 
    ]
    
    # Initialize structures
    fluxes = {}
    thermal = {}
    rad = {}

    for name in scalar_fields:
        if hasattr(fluxes, name) or name in ['Rntot', 'lEtot', 'Htot', 'Atot', 'Tcave', 'Tsave', 'raa', 'rawc', 'raws', 'ustar']:
            # Assign to the correct structure based on context/naming
            if name in ['Tcave', 'Tsave']:
                thermal[name] = np.nan
            elif name in ['raa', 'rawc', 'raws', 'ustar']:
                # These belong in the 'resistance' struct, but were put in 'thermal' here (as in MATLAB)
                thermal[name] = np.nan
            else:
                fluxes[name] = np.nan
        else:
            # Assume other scalars go to rad
            rad[name] = np.nan
    
    # Thermal array
    thermal['Ts'] = np.full(2, np.nan)

    # Spectral fields (size of wlF)
    wlF_size = spectral['wlF'].shape[0] if spectral['wlF'].ndim > 0 else 1
    rad['LoF_'] = np.full(wlF_size, np.nan).reshape(-1, 1)
    rad['Fhem_'] = np.full(wlF_size, np.nan).reshape(-1, 1)

    # Spectral fields (size of wlS)
    wlS_size = spectral['wlS'].shape[0] if spectral['wlS'].ndim > 0 else 1
    # Note: MATLAB does not correctly initialize Lo_ here, but we do for safety.
    rad['Lout_'] = np.full(wlS_size, np.nan).reshape(-1, 1)
    rad['Lo_'] = np.full(wlS_size, np.nan).reshape(-1, 1)
    
    thermal['Ta'] = np.nan
    
    return rad, thermal, fluxes