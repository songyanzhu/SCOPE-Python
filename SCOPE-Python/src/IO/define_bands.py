import numpy as np

def define_bands():
    """
    Translation of src/IO/define_bands.m
    Defines spectral regions for SCOPE.
    """

    spectral = {}

    # 3 spectral regions for SCOPE
    # MATLAB's ':' creates ranges. np.arange needs explicit step.
    reg1 = np.arange(400, 2401, 1)    # 400:1:2400
    reg2 = np.arange(2500, 15001, 100) # 2500:100:15000
    reg3 = np.arange(16000, 50001, 1000)# 16000:1000:50000

    # MATLAB: [reg1 reg2 reg3]' (column vector)
    spectral['wlS'] = np.concatenate((reg1, reg2, reg3)).reshape(-1, 1)

    # Other spectral (sub)regions
    spectral['wlP'] = reg1.reshape(-1, 1)                             # PROSPECT data range
    spectral['wlE'] = np.arange(400, 751, 1)                         # excitation in E-F matrix
    spectral['wlF'] = np.arange(640, 851, 1)                         # chlorophyll fluorescence in E-F matrix
    spectral['wlO'] = reg1                                           # optical part
    spectral['wlT'] = np.concatenate((reg2, reg3))                   # thermal part
    spectral['wlZ'] = np.arange(500, 601, 1)                         # xanthophyll region

    # PAR range
    wlS = spectral['wlS']
    spectral['wlPAR'] = wlS[(wlS >= 400) & (wlS <= 700)]

    # Indices (1-based in MATLAB, 0-based in Python. We use Python's convention, but the indices are calculated based on the MATLAB structure)
    spectral['IwlF'] = np.arange(640, 851, 1) - 400 # 640:1:850 nm, starting from 400nm (index 400 = 1st element)
    spectral['IwlP'] = np.arange(1, len(reg1) + 1)
    spectral['IwlT'] = np.arange(len(reg1) + 1, len(reg1) + len(reg2) + len(reg3) + 1)

    # Data used by aggreg routine to read in MODTRAN data
    spectral['SCOPEspec'] = {
        'nreg': 3,
        'start': np.array([400, 2500, 16000]),
        'end': np.array([2400, 15000, 50000]),
        'res': np.array([1, 100, 1000])
    }
    
    return spectral