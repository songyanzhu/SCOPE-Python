import numpy as np
import os
import math
from typing import Tuple, Dict

# Constants of Planck function (Wouter Verhoef Sept. 2011)
C1 = 1.191066e-22 # W * m^2 / sr / nm^-5
C2 = 14388.33    # K * nm

def soltir_tp7(filename: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Reads MODTRAN tp7 file and applies a new MIT algorithm to derive
    18 spectral functions for atmospheric correction and simulations at BOA and TOA.
    
    Args:
        filename: Path to the MODTRAN .tp7 output file.

    Returns:
        Tuple[wn, wl, Tall]: Wavenumber (cm^-1), Wavelength (nm), 
                             and the 18 derived atmospheric functions.
    """
    
    # 1. File and Header Parsing
    try:
        fid = open(filename, 'r')
    except FileNotFoundError:
        print(f"Error: File not found at {filename}")
        return np.array([]), np.array([]), np.array([])

    modname = os.path.splitext(filename)[0]
    outname = f"{modname}.atm"

    # Skip first 7 header lines
    for _ in range(7):
        fid.readline()

    # Read solar zenith angle (tts)
    rline = fid.readline().strip()
    # str2num is approximated by converting to float array
    rline_floats = [float(x) for x in rline.split() if x]
    if len(rline_floats) < 4:
        fid.close()
        raise ValueError("Could not read solar zenith angle (tts) from header.")
        
    tts = rline_floats[3]
    cts = np.cos(np.deg2rad(tts))

    # Read spectral range (wns, wne, wstep)
    s = fid.readline()
    # sscanf('%10f', 3) is approximated by simple splitting
    rline_floats = [float(s[i:i+10].strip()) for i in range(0, 30, 10)]
    wns, wne, wstep = rline_floats[0], rline_floats[1], rline_floats[2]
    
    # Calculate number of spectral points
    nw = int(np.round((wne - wns) / wstep)) + 1

    # Skip 2 lines before data start
    fid.readline()
    fid.readline()

    datarec = np.zeros((nw, 15, 4)) # MODTRAN5.2.1 tp7 15 column output format
    Tall = np.zeros((nw, 18))       # 18 output spectra

    # 2. Data Reading Loop (4 MODTRAN runs)
    for ipass in range(4):
        for il in range(nw):
            s = fid.readline()
            dline = [float(x) for x in s.split() if x]
            # Ensure 15 columns are read, padding with zeros if necessary (though MODTRAN standardizes output)
            datarec[il, :, ipass] = dline[:15]
        
        # Skip 12 lines of trailer info after each pass
        for j in range(12):
            fid.readline()

    fid.close()

    # 3. Core Calculations (MIT Algorithm)

    wn = datarec[:, 0, 0]
    fac = wn * wn
    wl = 1.0e7 / wn
    # wls = wl[nw - 1] # MATLAB's 1-based indexing for wl(nw)
    # wle = wl[0]

    # Transmittance
    tran_boa = datarec[:, 1, 0] # COL 2, Pass 1 (BOA, a=0.5)
    tran_toa = datarec[:, 1, 2] # COL 2, Pass 3 (TOA, a=0.0)
    too = tran_toa

    # TOA Sun Irradiance
    # datarec(:,14,4) is ESUN (TOA Incident Solar Flux)
    toasun = datarec[:, 13, 3] * fac / np.pi * cts  # = Eso cos(tts) / pi

    # --- BOA Derivatives (Pass 1 & 2) ---
    
    # Radiance reflected by ground (a=0.5)
    grfl50_boa = datarec[:, 6, 0] * fac # COL 7, Pass 1 (BOA, a=0.5)
    # Radiance reflected by ground (a=1.0)
    grfl100_boa = datarec[:, 6, 1] * fac # COL 7, Pass 2 (BOA, a=1.0)
    # Radiance reflected by surface emission (a=0.5)
    sfem50 = datarec[:, 3, 0] * fac # COL 4, Pass 1 (BOA, a=0.5)
    # Radiance of surface emission (a=0.0)
    sfem0_boa = 2 * sfem50
    
    # delgtot_boa = L(a=1) - L(a=0)
    delgtot_boa = grfl100_boa - sfem0_boa
    # crdd = L(a=0.5) - L_self_emission / (L(a=1) - L(a=0.5) - L_self_emission)
    crdd = (grfl50_boa - sfem50) / (grfl100_boa - grfl50_boa - sfem50)
    
    # rdd: directional-hemispherical reflectance factor for diffuse radiation at BOA (albedo)
    rdd = np.maximum(0, 1 - crdd)

    # OpT: Upwelling path radiance emitted by the ground
    OpT = crdd * delgtot_boa / tran_boa

    # OpT at 6500 nm for Planck reference (Nearest Neighbor Interpolation)
    wlp = 6500.0
    idx_wlp = np.argmin(np.abs(wl - wlp))
    Lp = OpT[idx_wlp]
    
    # Planck Brightness Temperature (Tp)
    # T = c2 / (lambda*1e-3 * log(1 + c1 * (lambda*1e-9)^(-5) / L))
    # wl is in nm, convert to m (1e-9) and cm-1 (1e-3)
    lambda_m = wlp * 1e-9
    lambda_cm_inv = wlp * 1e-3 
    
    Tp = C2 / (lambda_cm_inv * np.log(1 + C1 * (lambda_m)**(-5) / Lp))
    
    # Lmin: Minimum upwelling path radiance (used for correction)
    # Lmin = c1*(wl*1e-9).^(-5)./(exp(c2./(wl*1e-3*Tp))-1);
    wl_m = wl * 1e-9
    wl_cm_inv = wl * 1e-3
    Lmin = C1 * (wl_m)**(-5) / (np.exp(C2 / (wl_cm_inv * Tp)) - 1)


    # --- TOA Derivatives (Pass 3 & 4) ---
    
    # TOA Radiance reflected by ground (a=1.0)
    grfl100_toa = datarec[:, 6, 3] * fac # COL 7, Pass 4 (TOA, a=1.0)
    # TOA Radiance reflected by surface emission (a=0.0)
    sfem0 = datarec[:, 3, 2] * fac # COL 4, Pass 3 (TOA, a=0.0)

    # delgtot_toa = L(a=1) - L(a=0) at TOA
    delgtot_toa = grfl100_toa - sfem0
    # OpTtran = upwelling path radiance from atmosphere/surface emission 
    OpTtran = crdd * delgtot_toa
    
    # Path radiance (scattered)
    path100 = datarec[:, 4, 3] * fac # COL 5, Pass 4
    path0 = datarec[:, 4, 2] * fac # COL 5, Pass 3
    delpath = path100 - path0

    # Path radiance (thermal/atmospheric emission)
    ptem100 = datarec[:, 2, 3] * fac # COL 3, Pass 4
    ptem0 = datarec[:, 2, 2] * fac # COL 3, Pass 3
    delptem = np.maximum(0, ptem100 - ptem0)
    
    # delatmo = delpath + delptem
    delatmo = delpath + delptem

    # rso: TOA path radiance (Lp,0) / toasun (E_sun * cos(tts) / pi)
    rso = path0 / toasun

    # fO, fT: Fractions of scattered vs thermal emission in OpT (upwelling path radiance)
    iT = (wl > 4600) # Indices where thermal emission dominates (IR)
    # Indices where a correction is needed due to one component being zero
    # ia = (~iT & delpath == 0) | (iT & delptem == 0)
    ia = (~iT & (delpath == 0)) | (iT & (delptem == 0)) 
    
    # fO = delpath / delatmo * ~ia + ~iT * ia
    # Add a small epsilon to avoid division by zero where delatmo is zero
    epsilon = 1e-60
    fO = np.where(~ia, delpath / (delatmo + epsilon), np.where(~iT, 1.0, 0.0))
    # fT = delptem / delatmo * ~ia + iT * ia
    fT = np.where(~ia, delptem / (delatmo + epsilon), np.where(iT, 1.0, 0.0))
    
    # O & T: Upwelling path radiance due to scattering (O) and thermal emission (T)
    O = fO * OpT
    T = fT * OpT

    # Correct cases where T = 0 using Lmin
    i0 = (T == 0)
    T = np.where(i0, Lmin, T)
    fT = np.where(i0, T / OpT, fT)
    
    Ttran = fT * OpTtran
    Otran = OpTtran - Ttran

    # tdo: Downwelling diffuse transmittance (tdo)
    # tdo = delatmo / delgtot_toa * tran_toa
    tdo = delatmo / (delgtot_toa + 1e-6) * tran_toa
    tdo = np.where(ia, 0.0, tdo) # Force 0 where correction was applied

    # --- Direct Transmittance Components ---
    
    # gsun100_toa: Direct solar flux reaching TOA (Esun * cos(tts) * too)
    gsun100_toa = datarec[:, 7, 3] * fac # COL 8, Pass 4
    # gsun100_boa: Direct solar flux reaching BOA (Esun * cos(tts) * tss)
    gsun100_boa = datarec[:, 7, 1] * fac # COL 8, Pass 2

    # tsstoo: TOA direct transmittance (too)
    tsstoo = gsun100_toa / toasun
    # tss: BOA direct transmittance (tss)
    tss = gsun100_boa / toasun
    
    # tsdM: Downwelling diffuse transmittance (tdo) from the MODTRAN approximation
    # tsdM = O/toasun - tss
    tsdM = O / toasun - tss
    tsdM = np.where(ia, 0.0, tsdM)

    # --- Log Regression Prediction for tsd ---

    # Ir: Spectral regions for regression (Near IR windows)
    Ir = (wl > 1500) & (wl < 1700) | (wl > 2100) & (wl < 2300) | (wl > 3950) & (wl < 4100)
    
    y = np.log(tsdM[Ir] / tss[Ir]) - np.log(tdo[Ir] / too[Ir])
    x = np.log(wl[Ir])
    
    n = len(x)
    xm = np.sum(x) / n
    ym = np.sum(y) / n
    
    # Linear Regression: y = a*x + b
    # a = Cov(x, y) / Var(x)
    a = np.dot((x - xm), (y - ym)) / np.dot((x - xm), (x - xm))
    b = ym - a * xm
    
    # Prediction (p) and predicted diffuse transmittance (tsdp)
    p = a * np.log(wl) + b
    tsdp = tss * tdo / (too + 1e-60) * np.exp(p)

    # weight proportional to delatmo squared
    wgtM = delatmo**2

    # Combined tsd (weighted average)
    tsd = (wgtM * tsdM + tsdp) / (wgtM + 1)
    
    # fsun: fraction of direct/total solar radiation (direct / (direct + diffuse))
    fsun = (tss + 1e-60) / (tss + 1e-60 + tsd)
    
    # iH: Indices in the mid-IR where tsdtoo prediction is unreliable
    iH = (wl > 2600)

    # tsdtoo: TOA downwelling diffuse transmittance
    tsdtoo = Otran / toasun - tsstoo
    # Correction in the mid-IR (iH)
    tsdtoo = np.where(iH, tsd * too, tsdtoo)

    # --- The 18 Derived Functions (T1...T18) ---
    
    # T1: toasun
    # T2: rso (TOA atmospheric reflectance for diffuse)
    # T3: rdd (BOA hemispherical diffuse reflectance/albedo for direct)
    # T4: tss (BOA direct transmittance)
    # T5: tsd (BOA diffuse transmittance)
    # T6: too (TOA direct transmittance)
    # T7: tdo (TOA diffuse transmittance)
    
    tssrddtoo = tsstoo * rdd           # T8
    toordd = too * rdd                # T9
    tssrdd = tss * rdd                # T10
    tsstdo = fsun * delpath * crdd / toasun # T11
    tsdtdo = (1 - fsun) * delpath * crdd / toasun # T12
    
    # Lat: Upwelling thermal path radiance from atmosphere [mW m-2 sr-1 nm-1]
    Lat = ptem0 - sfem0_boa / crdd * tdo # T13
    
    # Lab: Upwelling surface emission + path radiance [mW m-2 sr-1 nm-1]
    Lab = T + sfem0_boa * crdd        # T14
    
    # Labtoo: Lab at TOA, direct component [mW m-2 sr-1 nm-1]
    Labtoo = Ttran + crdd * sfem0_boa * too # T15
    
    # Labtdo: Lab at TOA, diffuse component [mW m-2 sr-1 nm-1]
    Labtdo = crdd * (delptem + sfem0_boa * tdo) # T16
    
    # Remaining two components: 
    # T17: tsdtoo (TOA diffuse transmittance for diffuse)
    # T18: tsstoo (already defined, TOA direct transmittance)
    # Note: MATLAB code provides 18 columns, the labels suggest 18 parameters.
    # Rechecking the MATLAB concatenation:
    # [toasun rso rdd tss tsd too tdo tsstoo tsdtoo tsstdo tsdtdo tssrdd toordd tssrddtoo Lat Lab Labtoo Labtdo]
    # This is 1 + 17 = 18 components, so `tsdtoo` and `tsstoo` are included explicitly.
    
    Tall = np.column_stack([
        toasun, rso, rdd, tss, tsd, too, tdo, tsstoo, tsdtoo, tsstdo, tsdtdo, 
        tssrdd, toordd, tssrddtoo, Lat, Lab, Labtoo, Labtdo
    ])

    # 4. Write data to output file
    try:
        fid = open(outname, 'w')

        # Header lines
        str1=' WN (cm-1) WL (nm)  '
        str2='      T1       T2       T3        T4        T5        T6        T7      '
        str3='  T8        T9        T10       T11       T12       T13       T14      '
        str4=' T15          T16          T17          T18'
        str5='                       toasun     rso      rdd       tss       tsd       too         tdo      tsstoo    tsdtoo    tsstdo    tsdtdo    tssrdd      toordd  tssrddtoo     Lat          Lab         Labtoo       Labtdo'
        str_header = str1 + str2 + str3 + str4
        
        # Write headers
        fid.write(f"{str_header}\n")
        fid.write(f"{str5}\n\n")

        # Write data rows
        for i in range(nw):
            row = Tall[i, :]
            
            # WN and WL (9.2f, 10.3f)
            str_data = f"{wn[i]:9.2f}{wl[i]:10.3f}"
            
            # T1 (10.5f)
            str_data += f"{row[0]:10.5f}"
            
            # T2 - T14 (13 items, 10.6f)
            for j in range(1, 14):
                str_data += f"{row[j]:10.6f}"
            
            # T15 - T18 (4 items, 14.6e)
            for j in range(14, 18):
                str_data += f"{row[j]:14.6e}"
                
            fid.write(f"{str_data}\n")

    except Exception as e:
        print(f"Error writing output file {outname}: {e}")
        
    finally:
        # Close all files explicitly (like fclose('all') in MATLAB)
        if 'fid' in locals() and not fid.closed:
            fid.close()

    return wn, wl, Tall
