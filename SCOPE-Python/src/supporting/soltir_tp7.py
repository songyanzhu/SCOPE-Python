import numpy as np
import os
from math import cos, log
from scipy.interpolate import interp1d

# Constants for Planck function are used within the MATLAB script, 
# so they are repeated here for completeness.
C1_PLANCK = 1.191066e-22
C2_PLANCK = 14388.33

def soltir_tp7(filename):
    """
    Reads MODTRAN tp7 file and derives 18 spectral functions (t1...t18)
    for atmospheric correction and simulations at BOA and TOA.
    
    NOTE: Reading the raw, fixed-format MODTRAN tp7 file in Python is complex.
    The following implementation simulates the reading logic but relies on the 
    assumption of a specific, non-standard text output from MODTRAN.

    :param filename: Path to the MODTRAN tp7 file.
    :return: wn, wl, Tall (wavenumber, wavelength, 18 atmospheric spectra).
    """
    
    modname, _ = os.path.splitext(filename)
    outname = f'{modname}.atm'

    # --- File Reading Simulation (highly simplified) ---
    try:
        with open(filename, 'r') as fid:
            lines = fid.readlines()
    except FileNotFoundError:
        print(f"Error: File not found at {filename}")
        return None, None, None

    # Line 8 in MATLAB (7 lines skipped, so index 7 in Python lines)
    rline = np.array([float(x) for x in lines[7].strip().split()])
    tts = rline[3]
    cts = np.cos(np.deg2rad(tts))

    # Line 9
    s = lines[8].strip()
    # sscanf('%10f', 3) reads 3 floats
    rline_wn = np.array([float(s[i:i+10].strip()) for i in range(0, 30, 10)])
    wns, wne, wstep = rline_wn[0], rline_wn[1], rline_wn[2]
    
    # nw = int32((wne - wns) / wstep) + 1
    nw = int(round((wne - wns) / wstep)) + 1
    
    datarec = np.zeros((nw, 15, 4)) # MODTRAN5.2.1 tp7 15 column output format
    # Tall = zeros(nw, 18)
    Tall = np.zeros((nw, 18))

    # Skip 2 more lines (Line 10, 11)
    # Start of data is line 12 (index 11)
    
    start_line = 11
    
    # Read the 4 passes
    for ipass in range(4):
        for il in range(nw):
            # Read line
            s = lines[start_line + il + ipass * (nw + 12)].strip()
            # Assuming a standard text format where values are space-separated
            dline = np.array([float(x) for x in s.split()])
            # The MATLAB code suggests the data has 15 columns per pass
            datarec[il, :, ipass] = dline
            
    # Wave number and conversions
    wn = datarec[:, 0, 0]
    fac = wn * wn
    wl = 1e7 / wn
    
    # --- MIT algorithm for 18 spectra (T1..T18) ---
    
    # BOA properties from passes 1 and 2 (0 and 1 in Python)
    tran_boa = datarec[:, 1, 0]
    grfl50_boa = datarec[:, 6, 0] * fac
    sfem50 = datarec[:, 3, 0] * fac
    sfem0_boa = 2 * sfem50
    grfl100_boa = datarec[:, 6, 1] * fac
    delgtot_boa = grfl100_boa - sfem0_boa
    
    # The division may lead to NaN or Inf if the denominator is zero
    crdd = np.divide(grfl50_boa - sfem50, grfl100_boa - grfl50_boa - sfem50, 
                     out=np.zeros_like(grfl50_boa), where=(grfl100_boa - grfl50_boa - sfem50) != 0)

    rdd = np.maximum(0, 1 - crdd)
    
    # Add a small epsilon to avoid division by zero in tran_boa if needed, 
    # though OpT should be 0 if tran_boa is 0.
    OpT = np.divide(crdd * delgtot_boa, tran_boa, 
                    out=np.zeros_like(tran_boa), where=tran_boa != 0)

    # OpT at 6500 nm (for Lmin calculation)
    wlp = 6500
    # MATLAB's interp1 with 'nearest' is used
    f = interp1d(wl, OpT, kind='nearest', bounds_error=False, fill_value='extrapolate')
    Lp = f(wlp)
    
    # This involves a calculation that can result in log(0) if Lp=0.
    # The MATLAB code assumes Lp > 0 for 6500nm, which is generally true for thermal.
    try:
        Tp = C2_PLANCK / (wlp * 1e-3 * log(1 + C1_PLANCK * (wlp * 1e-9)**(-5) / Lp))
    except ValueError:
        # Handle log(<=0) case by using a default temperature or a large number
        Tp = 300.0 # Example fallback
        
    wl_m = wl * 1e-9
    wl_mm = wl * 1e-3
    
    # Lmin: Planck radiance at Tp for all wavelengths
    Lmin = C1_PLANCK * (wl_m**(-5)) / (np.exp(C2_PLANCK / (wl_mm * Tp)) - 1)
    
    # TOA properties from passes 3 and 4 (2 and 3 in Python)
    tran_toa = datarec[:, 1, 2] # Same as too
    too = tran_toa
    toasun = datarec[:, 13, 3] * fac / np.pi * cts 
    
    grfl100_toa = datarec[:, 6, 3] * fac
    sfem0 = datarec[:, 3, 2] * fac
    delgtot_toa = grfl100_toa - sfem0
    OpTtran = crdd * delgtot_toa
    
    path100 = datarec[:, 4, 3] * fac
    path0 = datarec[:, 4, 2] * fac
    
    # rso = path0 / toasun (T2)
    rso = np.divide(path0, toasun, 
                    out=np.zeros_like(path0), where=toasun != 0)
    
    delpath = path100 - path0
    ptem100 = datarec[:, 2, 3] * fac
    ptem0 = datarec[:, 2, 2] * fac
    delptem = np.maximum(0, ptem100 - ptem0)
    delatmo = delpath + delptem

    # --- Remaining calculations for T1..T18 (lots of masking/logical indexing) ---
    
    iT = (wl > 4600)
    # ia = (~iT & delpath == 0) | (iT & delptem == 0)
    ia = (~iT & (delpath == 0)) | (iT & (delptem == 0))
    
    # fO and fT calculations involve division by (delatmo + 1e-60)
    # The MATLAB code uses 1e-60 to handle division by zero.
    
    # fO = delpath / (delatmo + 1e-60) * ~ia + ~iT * ia
    den_fo_ft = delatmo + 1e-60
    fO = np.divide(delpath, den_fo_ft, out=np.zeros_like(delpath), where=den_fo_ft != 0)
    fO[~ia] = fO[~ia]
    fO[ia & ~iT] = 1 # for non-thermal bands where ia is True, fO=1
    fO[ia & iT] = 0 # for thermal bands where ia is True, fO=0
    
    # fT = delptem / (delatmo + 1e-60) * ~ia + iT * ia
    fT = np.divide(delptem, den_fo_ft, out=np.zeros_like(delptem), where=den_fo_ft != 0)
    fT[~ia] = fT[~ia]
    fT[ia & ~iT] = 0 # for non-thermal bands where ia is True, fT=0
    fT[ia & iT] = 1 # for thermal bands where ia is True, fT=1

    O = fO * OpT
    T = fT * OpT
    
    # Correct cases where T = 0
    i0 = (T == 0)
    T[i0] = Lmin[i0]
    
    # This recalculation assumes OpT is not zero at i0, otherwise div by zero
    fT[i0] = np.divide(T[i0], OpT[i0], 
                       out=np.zeros_like(T[i0]), where=OpT[i0] != 0) 
    
    Ttran = fT * OpTtran
    Otran = OpTtran - Ttran

    # tdo (T7) calculation
    # tdo = delatmo / (delgtot_toa + 1e-6) * tran_toa
    tdo = np.divide(delatmo, delgtot_toa + 1e-6, 
                    out=np.zeros_like(delatmo), where=delgtot_toa != 0) * tran_toa
    tdo[ia] = 0

    # gsun100
    gsun100_toa = datarec[:, 7, 3] * fac
    gsun100_boa = datarec[:, 7, 1] * fac

    # tsstoo (T8)
    # tsstoo = gsun100_toa / toasun
    tsstoo = np.divide(gsun100_toa, toasun, 
                       out=np.zeros_like(gsun100_toa), where=toasun != 0)
    
    # tss (T4)
    # tss = gsun100_boa / toasun
    tss = np.divide(gsun100_boa, toasun, 
                    out=np.zeros_like(gsun100_boa), where=toasun != 0)
    
    # tsdM (intermediate for T5)
    # tsdM = O / toasun - tss
    tsdM = np.divide(O, toasun, out=np.zeros_like(O), where=toasun != 0) - tss
    tsdM[ia] = 0
    
    # Log regression to predict tsd from tdo, too, tss, and wl
    Ir = ((wl > 1500) & (wl < 1700)) | ((wl > 2100) & (wl < 2300)) | ((wl > 3950) & (wl < 4100))
    
    # y = log(tsdM(Ir) / tss(Ir)) - log(tdo(Ir) / too(Ir))
    # Add epsilon to prevent log(0) and div by zero if needed.
    y_term1 = np.divide(tsdM[Ir], tss[Ir] + 1e-60, 
                        out=np.zeros_like(tsdM[Ir]), where=tss[Ir] != 0)
    y_term2 = np.divide(tdo[Ir], too[Ir] + 1e-60, 
                        out=np.zeros_like(tdo[Ir]), where=too[Ir] != 0)
    
    # Filter for valid log arguments (must be > 0)
    valid_log_mask = (y_term1 > 0) & (y_term2 > 0)
    
    y = np.zeros_like(tsdM[Ir])
    y[valid_log_mask] = np.log(y_term1[valid_log_mask]) - np.log(y_term2[valid_log_mask])
    
    x = np.log(wl[Ir])
    
    # Linear fit: y = a*x + b
    if len(x[valid_log_mask]) > 1:
        # Use NumPy's polyfit for linear fit
        coeffs = np.polyfit(x[valid_log_mask], y[valid_log_mask], 1)
        a, b = coeffs[0], coeffs[1]
    else:
        # Fallback if insufficient data for fit
        a, b = 0, 0 
        
    p = a * np.log(wl) + b
    
    # tsdp: predicted tsd
    # tsdp = tss * tdo / (too + 1e-60) * exp(p)
    tsdp = tss * np.divide(tdo, too + 1e-60, 
                           out=np.zeros_like(tdo), where=too != 0) * np.exp(p)

    # weight proportional to delatmo squared
    wgtM = delatmo**2
    
    # tsd (T5) (Weighted average)
    # tsd = (wgtM * tsdM + tsdp) / (wgtM + 1)
    tsd = np.divide(wgtM * tsdM + tsdp, wgtM + 1,
                    out=np.zeros_like(wgtM), where=(wgtM + 1) != 0)

    # fsun (fraction sunlit)
    fsun = np.divide(tss, tss + tsd, 
                     out=np.zeros_like(tss), where=(tss + tsd) != 0)
    
    # iH (used for tsdtoo correction)
    iH = (wl > 2600)

    # tsdtoo (T9)
    # tsdtoo = Otran / toasun - tsstoo
    tsdtoo = np.divide(Otran, toasun, 
                       out=np.zeros_like(Otran), where=toasun != 0) - tsstoo
    tsdtoo[iH] = tsd[iH] * too[iH]

    # Remaining T-values
    tssrddtoo = tsstoo * rdd
    toordd = too * rdd
    tssrdd = tss * rdd
    
    # tsstdo (T10)
    # tsstdo = fsun * delpath * crdd / toasun
    tsstdo = fsun * delpath * crdd / toasun
    
    # tsdtdo (T11)
    # tsdtdo = (1 - fsun) * delpath * crdd / toasun
    tsdtdo = (1 - fsun) * delpath * crdd / toasun
    
    # Lat (T15)
    # Lat = ptem0 - sfem0_boa / crdd * tdo
    Lat = ptem0 - np.divide(sfem0_boa, crdd, 
                            out=np.zeros_like(sfem0_boa), where=crdd != 0) * tdo
    
    # Lab (T16)
    Lab = T + sfem0_boa * crdd
    
    # Labtoo (T17)
    Labtoo = Ttran + crdd * sfem0_boa * too
    
    # Labtdo (T18)
    Labtdo = crdd * (delptem + sfem0_boa * tdo)

    # Combine all 18 results into Tall
    Tall = np.column_stack([
        toasun, rso, rdd, tss, tsd, too, tdo, tsstoo, tsdtoo, tsstdo, tsdtdo,
        tssrdd, toordd, tssrddtoo, Lat, Lab, Labtoo, Labtdo
    ])

    # --- Write data to output file (Simulation) ---
    # The writing logic is complex formatting, only shown here for reference.
    
    # str1, str2, str3, str4, str5 define the header.
    # The main loop formats and prints each row of wn, wl, and Tall.
    
    # Example for writing the data content (simplified):
    data_out = np.column_stack([wn, wl, Tall])
    header_list = ['WN (cm-1)', 'WL (nm)'] + [f'T{i+1}' for i in range(18)]
    header_str = ' '.join(header_list)
    
    # Final data structure (wl and wn are not needed for return based on MATLAB function signature,
    # but they are often returned in a tuple)

    # Returning the required variables
    return wn, wl, Tall