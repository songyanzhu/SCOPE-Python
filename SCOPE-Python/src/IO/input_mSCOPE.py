import numpy as np
import pandas as pd
import warnings
import os

def input_mSCOPE(parameter_file):
    """
    Translation of src/IO/input_mSCOPE.m
    Loads multilayer input parameters from a CSV file (analogous to MATLAB's xlsread).
    """
    try:
        # MATLAB: mlayer = xlsread(parameter_file)
        # Assuming the CSV equivalent contains numbers only
        df = pd.read_csv(parameter_file, header=None).fillna(0)
        mlayer = df.values # Get NumPy array
    except Exception as e:
        raise IOError(f"Failed to read mSCOPE parameter file {parameter_file}: {e}")

    mly = {}
    
    # MATLAB indexing is 1-based, we use 0-based
    # mlayer(1,1) is the number of layers
    mly['nly'] = int(mlayer[0, 0])
    
    # mlayer(row, :) gives the row as a vector (MATLAB)
    # The first column is often a descriptor, so the layer values start from the second column (index 1)
    # The MATLAB implementation seems to read the entire row and then sum it later,
    # implicitly assuming the parameters start from column 1 or 2. We assume the MATLAB 
    # indexing corresponds to rows in the sheet, and the data is column-wise.

    # Assuming parameters start from row 2 (index 1 in 0-based) for LAI:
    mly['pLAI'] = mlayer[1, :]
    mly['pCab'] = mlayer[2, :]
    mly['pCca'] = mlayer[3, :]
    mly['pCdm'] = mlayer[4, :]
    mly['pCw'] = mlayer[5, :]
    mly['pCs'] = mlayer[6, :]
    mly['pN'] = mlayer[7, :]
    
    # Slice to keep only the actual layers if the rows are longer
    if mly['pLAI'].size > mly['nly']:
        mly['pLAI'] = mly['pLAI'][:mly['nly']]
        mly['pCab'] = mly['pCab'][:mly['nly']]
        mly['pCca'] = mly['pCca'][:mly['nly']]
        mly['pCdm'] = mly['pCdm'][:mly['nly']]
        mly['pCw'] = mly['pCw'][:mly['nly']]
        mly['pCs'] = mly['pCs'][:mly['nly']]
        mly['pN'] = mly['pN'][:mly['nly']]
        
    # totLAI = sum(mly.pLAI)
    mly['totLAI'] = np.sum(mly['pLAI'])
    
    # MATLAB only requires the total LAI and layer profiles (pLAI, pCab, etc.)
    return mly