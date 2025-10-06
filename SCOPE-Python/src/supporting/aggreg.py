import numpy as np
import pandas as pd
from scipy.io import loadmat # MATLAB's importdata sometimes handles .mat files

def aggreg(atmfile, SCOPEspec):
    """
    Aggregate MODTRAN data over SCOPE bands by averaging (over rectangular band passes).

    :param atmfile: Path to the MODTRAN .atm file (assumed to be text/CSV data).
    :param SCOPEspec: Structure/dict with nreg, start, end, res.
    :return: Aggregated data matrix M.
    """
    
    # MATLAB's `importdata` is complex; assuming the .atm file is a simple
    # space-delimited text file (typical for MODTRAN output) where 
    # the second column is wlM and columns 3-20 are T.
    try:
        # Use pandas to read the data, skipping potential header rows
        # Assuming header has been skipped by MATLAB's importdata
        s = pd.read_csv(atmfile, sep='\s+', header=None, skiprows=4) # Example skip, adjust as needed
        # Assuming MATLAB's `s.data` corresponds to the dataframe's numerical content
        data = s.values
    except Exception as e:
        print(f"Error reading {atmfile}: {e}. Ensure file format is correct.")
        return None

    # wlM = s.data(:,2); 
    # MATLAB uses 1-based indexing and the 2nd column. Python uses 0-based.
    wlM = data[:, 1]
    
    # T = s.data(:,3:20);
    T = data[:, 2:20]

    # Extract 6 relevant columns from T (MATLAB 1-based index to Python 0-based)
    # T column indices (MATLAB): 1, 3, 4, 5, 12, 16
    # Corresponding Python 0-based indices: 0, 2, 3, 4, 11, 15
    U = T[:, [0, 2, 3, 4, 11, 15]]

    nwM = len(wlM)

    nreg = SCOPEspec['nreg']
    streg = SCOPEspec['start']
    enreg = SCOPEspec['end']
    width = SCOPEspec['res']

    # Convert to NumPy arrays for consistent slicing
    streg = np.asarray(streg)
    enreg = np.asarray(enreg)
    width = np.asarray(width)
    
    # Nr. of bands in each region
    # MATLAB's int32 division rounds towards zero, Python's `//` does floor division.
    # The MATLAB expression is: int32((enreg-streg)./width) + 1
    nwreg = np.floor((enreg - streg) / width).astype(int) + 1

    # off = int32(zeros(nreg,1));
    off = np.zeros(nreg, dtype=int)

    for i in range(1, nreg):
        off[i] = off[i-1] + nwreg[i-1]

    nwS = np.sum(nwreg)
    n = np.zeros(nwS, dtype=int)  # Count of MODTRAN data contributing to a band
    S = np.zeros((nwS, 6))      # Initialize sums (float)

    # j = int32(zeros(nreg,1)); # Band index within regions
    j = np.zeros(nreg, dtype=int) 

    for iwl in range(nwM):
        w = wlM[iwl]  # MODTRAN wavelength
        for r in range(nreg):
            # j(r) = int32(round(w-streg(r))./(width(r)))+1;
            # MATLAB round is to nearest integer. Python's `round` for a single number is similar.
            # Then +1 for 1-based indexing.
            # Python's round on a single value:
            j_r = int(round((w - streg[r]) / width[r])) + 1 
            j[r] = j_r
            
            # test if index is in valid range
            if j[r] > 0 and j[r] <= nwreg[r]:
                k = j[r] + off[r] - 1 # SCOPE band index (adjusting for Python 0-based k)
                
                # S(k,:) = S(k,:) + U(iwl,:); 
                S[k, :] = S[k, :] + U[iwl, :]
                
                # n(k) = n(k) + 1;
                n[k] = n[k] + 1

    # M = zeros(size(S,1),6);
    M = np.zeros(S.shape)
    
    # Calculate averages per SCOPE band
    # Handle division by zero for bands with no contributing data
    # np.divide handles division by zero by inserting Inf/NaN, then we can use np.where
    M = np.divide(S, n[:, np.newaxis], out=np.zeros_like(S, dtype=float), where=n[:, np.newaxis] != 0)

    # A simpler way, assuming bands with n=0 should result in an average of 0.
    # for i in range(6):
    #     M[:, i] = np.where(n != 0, S[:, i] / n, 0)
        
    return M