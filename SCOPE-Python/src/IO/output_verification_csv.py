import numpy as np
import os
import pandas as pd
import warnings
# Matplotlib is the standard Python plotting library.
try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

# --- Helper Function (Load CSV) ---
def _dlmread_mock(file_path, delimiter=',', header_rows=2):
    """
    Mocks MATLAB's dlmread(..., ',', 2, 0) or dlmread(..., ',', 1, 0)
    Reads CSV, skipping initial header/unit rows.
    """
    try:
        # Python's pandas can handle this skip
        df = pd.read_csv(file_path, delimiter=delimiter, header=None, skiprows=header_rows)
        return df.values
    except Exception as e:
        warnings.warn(f"Failed to read data from {file_path}: {e}")
        return None

# --- Main Function ---

def output_verification_csv(Output_dir, verification_dir):
    """
    Translation of src/IO/output_verification_csv.m
    Compares the latest output CSV files against a verification dataset.
    """

    # --- Setup Paths ---
    path0 = os.path.join('output', verification_dir)
    path1 = Output_dir

    if not os.path.isdir(path0):
        warnings.warn(f"Verification directory not found: {path0}")
        return

    # Get list of CSV files
    info0 = [f for f in os.listdir(path0) if f.endswith('.csv')]
    info1 = [f for f in os.listdir(path1) if f.endswith('.csv')]
    
    differentsize, differentcontent, differentnumberoffiles = False, False, False

    # --- Check Number of Files ---
    if len(info0) != len(info1):
        print(f"\nWarning: in the output file, {len(info1)} files were stored,")
        print(f"whereas there should be {len(info0)} files in this directory.")
        print("Check the simulation options that are specified the options tab of the input spreadsheet.")
        differentnumberoffiles = True

    # --- Check File Size and Content ---
    for n0 in info0:
        p0 = os.path.join(path0, n0)
        
        try:
            s0 = os.path.getsize(p0)
        except OSError:
            warnings.warn(f"Could not get size for verification file {p0}")
            continue

        if n0 in info1:
            p1 = os.path.join(path1, n0)
            try:
                s1 = os.path.getsize(p1)
            except OSError:
                warnings.warn(f"Could not get size for output file {p1}")
                continue

            # Check size
            if s0 != s1:
                print(f"\nWarning: the file size of {n0} is different from the verification output")
                print(f"({s1} instead of {s0} bytes)")
                differentsize = True
            
            # Check content
            header_rows = 2 # Default: skip header and units
            if n0 == 'pars_and_input_short.csv':
                continue # Skip short parameter file content check
            elif n0 == 'pars_and_input.csv':
                header_rows = 1 # Only skip header
                
            D0 = _dlmread_mock(p0, header_rows=header_rows)
            D1 = _dlmread_mock(p1, header_rows=header_rows)
            
            if D0 is not None and D1 is not None:
                if D0.shape == D1.shape:
                    # Sum of squared differences check (MATLAB: sum(sum(D0-D1).^2))>1E-9)
                    if np.sum(np.square(D0 - D1)) > 1E-9:
                        print(f"\nWarning: data in the output file {n0} are different from the verification output")
                        differentcontent = True
                        
                        # Plotting differences (if matplotlib is available)
                        if plt:
                            # Plotting logic is complex; we'll plot the first 4 variables as an example
                            fig_n = D0.shape[1] # Number of columns
                            spn = int(np.ceil(np.sqrt(min(fig_n, 9)))) # Limit to 9 plots
                            
                            plt.figure()
                            plt.suptitle(f"Differences in {n0}", fontsize=16)
                            
                            for m in range(min(fig_n, 9)):
                                plt.subplot(spn, spn, m + 1)
                                plt.plot(D0[:, m], 'k', label='Verification')
                                plt.plot(D1[:, m], 'r', label='Latest Run')
                                plt.title(f"Column {m+1}")
                            plt.show()

                else:
                    warnings.warn(f"File {n0} shapes do not match: {D0.shape} vs {D1.shape}")
            
        else:
            print(f"\nWarning: the file {n0} was not found in the output")
            
    # --- Final Summary ---
    if differentsize:
        print('\nWarning The size of some of the output files is different from the verification data')
        if differentcontent:
            print('Check if the startdate and enddate in the spreadsheet')
            print('and the verification data in are specified in "Dataset_Dir" in the Filenames tab of the input data spreadsheet')
        else:
            print('but the numerical values are the same. Possible cause: Different Matlab version')
            
    if not (differentsize or differentcontent or differentnumberoffiles):
        print('The output is the same as in the verification data set')