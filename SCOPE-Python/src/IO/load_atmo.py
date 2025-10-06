import numpy as np
import os

def load_atmo(atmfile, SCOPEspec):
    """
    Translation of src/IO/load_atmo.m
    Loads atmospheric input (irradiance spectra).
    """
    
    if not os.path.exists(atmfile):
        raise FileNotFoundError(f"Atmospheric file `{atmfile}` does not exist")

    # MATLAB: [~, ~, ext] = fileparts(atmfile);
    _, ext = os.path.splitext(atmfile)
    ext = ext.lower()

    atmo = {}

    if ext == '.atm':
        # MATLAB: atmo.M = aggreg(atmfile, SCOPEspec);
        # NOTE: The 'aggreg' function is not provided and is complex (MODTRAN data aggregation).
        # This will be replaced by a placeholder or raises an error in a functional setup.
        raise NotImplementedError("Reading of .atm files requires the 'aggreg' function, which is not implemented.")
    else:
        # Assume simple two-column CSV/TXT file: [Esun_ Esky_]
        try:
            raddata = np.loadtxt(atmfile)
        except Exception as e:
            raise IOError(f"Failed to load atmospheric data from {atmfile}: {e}")
            
        # MATLAB: atmo.Esun_ = raddata(:,1); atmo.Esky_ = raddata(:,2);
        if raddata.ndim == 1:
            # Handle single row case
            atmo['Esun_'] = raddata[0]
            atmo['Esky_'] = raddata[1]
        else:
            atmo['Esun_'] = raddata[:, 0]
            atmo['Esky_'] = raddata[:, 1]
            
    return atmo