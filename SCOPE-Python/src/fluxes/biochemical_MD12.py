import numpy as np

# --- Helper Functions (Private to Module) ---

def _gs_fun(Cs, RH, A, BallBerrySlope, BallBerry0):
    """Helper function for BallBerry to compute stomatal conductance gs (replicated from biochemical.py)."""
    gs = np.maximum(BallBerry0, BallBerrySlope * A * RH / (Cs + 1e-9) + BallBerry0)
    gs = np.where(np.isnan(Cs), np.nan, gs)
    return gs

def _ball_berry(Cs, RH, A, BallBerrySlope, BallBerry0, minCi, Ci_input=None):
    """
    MATLAB's BallBerry Model translation (non-iterative only in this file's context).
    """
    if Ci_input is not None and Ci_input.size > 0:
        Ci = Ci_input
        gs = None
        if A is not None and A.size > 0 and A.ndim == Ci_input.ndim:
            gs = _gs_fun(Cs, RH, A, BallBerrySlope, BallBerry0)
    elif np.all(BallBerry0 == 0) or A is None or A.size == 0 or A.ndim == 0:
        Ci = np.maximum(minCi * Cs, Cs * (1 - 1.6 / (BallBerrySlope * RH)))
        gs = None
    else:
        gs = _gs_fun(Cs, RH, A, BallBerrySlope, BallBerry0)
        Ci = np.maximum(minCi * Cs, Cs - 1.6 * A / gs)
    return Ci, gs

def satvap(T_c):
    """Saturated vapour pressure (hPa) at temperature T (Celsius)"""
    return 6.107 * 10**(7.5 * T_c / (237.3 + T_c))

def _md12(ps, Ja, Jms, kps, kf, kds, kDs):
    """
    MD12 algorithm for the computation of fluorescence yield.
    """
    # fs1: CO2-limited
    fs1 = ps * (kf / kps) / (1 - Ja / Jms)
    
    # fs2: light-limited
    par1 = kps / (kps - kds)                            
    par2 = par1 * (kf + kDs + kds) / kf                  
    fs2 = (par1 - ps) / par2                          

    fs = np.minimum(fs1, fs2)                              
    return fs

# --- Main Function ---

def biochemical_MD12(leafbio, meteo, options, constants, fV, Q_input=None):
    """
    Translation of src/fluxes/biochemical_MD12.m
    """
    # --- Input Parsing and Setup ---

    p = meteo['p'] * 1e2 # [hPa] -> [Pa]
    BallBerrySlope, BallBerry0 = leafbio['BallBerrySlope'], leafbio['BallBerry0']
    O, Cs, T, eb = meteo['Oa'], meteo['Cs'], meteo['T'], meteo['eb']
    Type, Tyear, beta, qLs, NPQs = leafbio['Type'], leafbio['Tyear'], leafbio['beta'], leafbio['qLs'], leafbio['kNPQs']
    Vcmax25, Rdparam = fV * leafbio['Vcmax25'], leafbio['RdPerVcmax25']
    
    Q = Q_input if Q_input is not None else meteo['Q']
    Q = np.where(Q == 0, 1e-9, Q)

    R, rhoa, Mair = constants['R'], constants['rhoa'], constants['Mair']

    # --- Unit Conversion ---
    T_celsius = T
    T = T + 273.15 * (T < 100) # [K]
    RH = np.minimum(1, eb / satvap(T_celsius))
    
    # Cs, O [bar] (p * 1e-11 for Cs, p * 1e-8 for O)
    Cs = Cs * p * 1e-11 
    O = O * p * 1e-8 
    
    # --- Photosynthetic Parameters @ TREF ---
    TREF, SCOOP = 25 + 273.15, 2862.
    Rdopt = Rdparam * Vcmax25
    
    if Type.lower() == 'c3':
        Jmo = Vcmax25 * 2.68
    else: # C4
        Jmo = Vcmax25 * 40/6 
        Vpmo, Vpr, x, alpha = Vcmax25 * 2.33, 80, 0.4, 0
        gbs = (0.0207 * Vcmax25 + 0.4806) * 1000. 

    # --- Temperature Corrections ---
    HARD, HAGSTAR = 46.39 * 1000, 37.83 * 1000
    CRD, CGSTAR = HARD / (R * TREF), HAGSTAR / (R * TREF)
    dum1, dum2 = R * T / 1000, R * TREF / 1000
    TDP = leafbio['TDP']
    
    Rd = Rdopt * np.exp(CRD - HARD / dum1)
    SCO = SCOOP / np.exp(CGSTAR - HAGSTAR / dum1)
    
    # Jmax, Vcmax (using Kattge and Knorr/Massad temp responses)
    HAJ, HDJ, DELTASJ = TDP['HAJ'] * 1000, TDP['HDJ'] * 1000, TDP['DELTASJ'] * 1000
    HAVCM, HDVC, DELTASVC = TDP['HAVCM'] * 1000, TDP['HDVC'] * 1000, TDP['DELTASVC'] * 1000
    
    Jmax = Jmo * np.exp(HAJ * (T - TREF) / (TREF * dum1)) * (1. + np.exp((TREF * DELTASJ - HDJ) / dum2)) / (1. + np.exp((T * DELTASJ - HDJ) / dum1))
    Vcmax = Vcmax25 * np.exp(HAVCM * (T - TREF) / (TREF * dum1)) * (1 + np.exp((TREF * DELTASVC - HDVC) / dum2)) / (1 + np.exp((T * DELTASVC - HDVC) / dum1))

    if Type.lower() == 'c3':
        HAKC, HAKO = TDP['HAKC'] * 1000, TDP['HAKO'] * 1000
        KCOP, KOOP = 404.9, 278.4
        CKC, CKO = HAKC / (R * TREF), HAKO / (R * TREF)
        Kc = KCOP * np.exp(CKC - HAKC / dum1) * 1e-11 * p 
        Ko = KOOP * np.exp(CKO - HAKO / dum1) * 1e-8 * p 
    else: # C4
        HAVPM, HDVP, DELTASVP = TDP['HAVPM'] * 1000, TDP['HDVP'] * 1000, TDP['DELTASVP'] * 1000
        KCOP, KOOP, KPOP = 944., 633., 82.
        Q10KC, Q10KO, Q10KP = 2.1, 1.2, 2.1
        
        Vpmax = Vpmo * np.exp(HAVPM * (T - TREF) / (TREF * dum1)) * (1 + np.exp((TREF * DELTASVP - HDVP) / dum2)) / (1 + np.exp((T * DELTASVP - HDVP) / dum1))
        Kc = KCOP * np.power(Q10KC, (T - TREF) / 10.0) * 1e-11 * p
        Ko = KOOP * np.power(Q10KO, (T - TREF) / 10.0) * 1e-8 * p
        Kp = KPOP * np.power(Q10KP, (T - TREF) / 10.0) * 1e-11 * p

    # --- Electron Transport and Fluorescence Parameters ---
    kf, kD, kd, po0max = 3.e7, 1.e8, 1.95e8, 0.88
    kPSII = (kD + kf) * po0max / (1. - po0max) 
    fo0 = kf / (kf + kPSII + kD)
    
    kps, kNPQs, kds = kPSII * qLs, NPQs * (kf + kD), kd * qLs
    kDs = kD + kNPQs
    Jms = Jmax * qLs
    po0 = kps / (kps + kf + kDs)
    THETA = (kps - kds) / (kps + kf + kDs)

    # J: Electron transport rate
    Q2 = beta * Q * po0
    J = (Q2 + Jms - np.sqrt((Q2 + Jms)**2 - 4 * THETA * Q2 * Jms)) / (2 * THETA)

    # --- Net Photosynthesis (A) ---
    minCi = 0.3 if Type.lower() == 'c3' else 0.1
    
    # Ci from non-iterative Ball-Berry (A=0 implicit assumption)
    Ci, _ = _ball_berry(Cs, RH, None, BallBerrySlope, 0, minCi) 

    if Type.lower() == 'c3':
        GSTAR, Cc = 0.5 * O / SCO, Ci
        Wc = Vcmax * Cc / (Cc + Kc * (1 + O / Ko))
        Wj = J * Cc / (4.5 * Cc + 10.5 * GSTAR)
        W = np.minimum(Wc, Wj)
        Ag = (1 - GSTAR / Cc) * W
        A = Ag - Rd
        
        # Ja: Actual linear electron transport rate (Wj > 0 check)
        Ja = np.where(Wj > 0, J * W / Wj, 0.0)
        
    else: # C4
        Cm, Rs, Rm, gam = Ci, 0.5 * Rd, 0.5 * Rd, 0.5 / SCO
        Vpc = Vpmax * Cm / (Cm + Kp)
        Vp = np.minimum(Vpc, Vpr)
        
        # Ac (CO2-limited) via quadratic
        dum1, dum2 = alpha / 0.047, Kc / Ko
        dum3, dum4, dum5, dum6 = Vp - Rm + gbs * Cm, Vcmax - Rd, gbs * Kc * (1 + O / Ko), gam * Vcmax
        a_ac, b_ac = 1.0 - dum1 * dum2, -(dum3 + dum4 + dum5 + dum1 * (dum6 + Rd * dum2))
        c_ac = dum4 * dum3 - dum6 * gbs * O + Rd * dum5
        Ac = (-b_ac - np.sqrt(b_ac**2 - 4 * a_ac * c_ac)) / (2.0 * a_ac)
        
        # Aj (Light-limited) via quadratic
        dum7, dum8 = x * J / 2.0 - Rm + gbs * Cm, (1.0 - x) * J / 3.0
        dum9, dum10 = dum8 - Rd, dum8 + Rd * 7.0 / 3.0
        a_aj, b_aj = 1.0 - 7.0 / 3.0 * gam * dum1, -(dum7 + dum9 + gbs * gam * O * 7.0 / 3.0 + dum1 * gam * dum10)
        c_aj = dum7 * dum9 - gbs * gam * O * dum10
        Aj = (-b_aj - np.sqrt(b_aj**2 - 4 * a_aj * c_aj)) / (2.0 * a_aj)
        
        A = np.minimum(Ac, Aj)
        Ja = J
        
        # Recalculate Ja if CO2-limited (A=Ac)
        ind_co2_limited = (A == Ac)
        if np.any(ind_co2_limited):
            A_ind = A[ind_co2_limited]
            A_safe = np.where(A_ind == 0, 1e-9, A_ind)
            
            # Complex quadratic for Ja
            a_ja = x * (1.0 - x) / 6.0 / A_safe
            b_ja = ((1.0 - x) / 3.0 * (gbs[ind_co2_limited] / A_safe * (Cm[ind_co2_limited] - Rm[ind_co2_limited] / gbs[ind_co2_limited] - gam * O[ind_co2_limited]) - 1.0 - alpha * gam / 0.047) 
                    - x / 2.0 * (1.0 + Rd[ind_co2_limited] / A_safe))
            c_ja = ((1.0 + Rd[ind_co2_limited] / A_safe) * (Rm[ind_co2_limited] - gbs[ind_co2_limited] * Cm[ind_co2_limited] - 7.0 * gbs[ind_co2_limited] * gam * O[ind_co2_limited] / 3.0) 
                    + (Rd[ind_co2_limited] + A_safe) * (1.0 - 7.0 * alpha * gam / 3.0 / 0.047))
            
            Ja_new = (-b_ja + np.sqrt(b_ja**2 - 4 * a_ja * c_ja)) / (2.0 * a_ja)
            Ja = np.where(ind_co2_limited, Ja_new, Ja)

    # --- Fluorescence and Resistances ---
    ps = Ja / (beta * Q)
    fs = _md12(ps, Ja, Jms, kps, kf, kds, kDs)
    eta = fs / fo0

    ppm2bar = 1e-6 * (p * 1e-5) # used in gs calculation
    gs = 1.6 * A * ppm2bar / (Cs - Ci)
    rcw = (rhoa / (Mair * 1e-3)) / gs
    rcw = np.where((A <= 0) & (rcw != 0), 0.625 * 1e6, rcw)

    # Ci back to ppm
    Ci_ppm = Ci / ppm2bar

    # --- Collect outputs ---
    biochem_out = {
        'A': A, 'Ci': Ci_ppm, 'ps': ps, 'eta': eta, 'fs': fs, 'rcw': rcw,
        'qE': np.full_like(rcw, np.nan),
        'Kn': NPQs + 0 * rcw,
        'Phi_N': kNPQs / (kNPQs + kD + kf + kps) + 0 * rcw,
        'Ja': Ja
    }

    return biochem_out