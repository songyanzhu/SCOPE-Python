import numpy as np

# --- Helper Functions (Private to Module) ---

def _sel_root(a, b, c, dsign):
    """
    MATLAB's sel_root - select a root based on dsign for the eqn ax^2 + bx + c = 0.
    """
    # Handle scalar a=0 case
    if np.isscalar(a) and a == 0:
        x = -c / b
    else:
        # dsign=0 is treated as dsign=-1 for selecting the smaller root
        dsign = np.where(dsign == 0, -1, dsign)
        # Calculate the discriminant: b^2 - 4ac
        discriminant = np.sqrt(b**2 - 4 * a * c)
        # Quadratic formula: (-b +/- sqrt(disc)) / 2a
        x = (-b + dsign * discriminant) / (2 * a)
    return x

def _gs_fun(Cs, RH, A, BallBerrySlope, BallBerry0):
    """Helper function to compute stomatal conductance gs."""
    # gs = max(BallBerry0, BallBerrySlope * A * RH / (Cs + 1e-9) + BallBerry0)
    gs = np.maximum(BallBerry0, BallBerrySlope * A * RH / (Cs + 1e-9) + BallBerry0)
    # Correct for NaN if Cs was NaN (MATLAB behavior)
    gs = np.where(np.isnan(Cs), np.nan, gs)
    return gs

def _ball_berry(Cs, RH, A, BallBerrySlope, BallBerry0, minCi, Ci_input=None):
    """
    MATLAB's BallBerry Model translation.
    """
    if Ci_input is not None and Ci_input.size > 0:
        # Ci is given (used in iterative solver: _ci_next -> _ball_berry(..., Ci_input))
        Ci = Ci_input
        gs = None
        if A is not None and A.size > 0 and A.ndim == Ci_input.ndim:
            gs = _gs_fun(Cs, RH, A, BallBerrySlope, BallBerry0)
    elif np.all(BallBerry0 == 0) or A is None or A.size == 0 or A.ndim == 0:
        # Non-iterative case (BallBerry0 = 0 or initial guess)
        # Ci = Cs * ( 1 - 1.6 / (m * RH) )
        Ci = np.maximum(minCi * Cs, Cs * (1 - 1.6 / (BallBerrySlope * RH)))
        gs = None
    else:
        # Iterative case (BallBerry0 > 0 and A is known from previous Ci guess)
        gs = _gs_fun(Cs, RH, A, BallBerrySlope, BallBerry0)
        # Ci = Cs - 1.6 * A/gs
        Ci = np.maximum(minCi * Cs, Cs - 1.6 * A / gs)
    
    return Ci, gs

def _fixedp_brent_ari(fun, x0, corner, tol):
    """
    MOCK IMPLEMENTATION of fixedp_brent_ari (vectorized fixed-point iteration).
    This mock uses a simple fixed-point iteration loop, which is less robust 
    than the original, but serves as a functional placeholder for the structure.
    """
    Ci_in = x0
    
    # Simple, vector-compatible fixed-point iteration (Brent's method is more complex)
    for _ in range(100): # Max 100 iterations (same as ebal's max counter)
        # fun returns (error, Ci_next_val). We need the Ci_next_val.
        err, Ci_next_val = fun(Ci_in)
        
        # Check for convergence
        if np.max(np.abs(err)) < tol:
            return Ci_next_val # Converged
            
        # Update guess
        Ci_in = Ci_next_val 
        
        # Simple damping for stability if not converging (not strictly in fixedp_brent_ari)
        if _ > 20:
             Ci_in = 0.9 * Ci_in + 0.1 * x0 
             
    # Return best effort if it didn't converge
    return Ci_next_val 

def _ci_next(Ci_in, Cs, RH, minCi, BallBerrySlope, BallBerry0, A_fun, ppm2bar):
    """
    Test-function for iteration: computes the difference (err) and the new estimate (Ci_out).
    """
    # Compute A using the guessed Ci_in
    A, _ = A_fun(Ci_in)
    A_bar = A * ppm2bar
    
    # Ci_out is the new estimate of Ci from the Ball-Berry equation (Ci = f(A))
    Ci_out, _ = _ball_berry(Cs, RH, A_bar, BallBerrySlope, BallBerry0, minCi) 

    # Error function for fixed-point iteration: err = f(x) - x
    err = Ci_out - Ci_in
    return err, Ci_out

# Use a global list/dictionary to mimic MATLAB's 'persistent' variable scope.
_compute_a_fcount = [0]

def _compute_a(Ci, Type, g_m, Vs_C3, MM_consts, Rd, Vcmax, Gamma_star, Je, effcon, atheta, kpepcase):
    """
    Compute Assimilation. Mimics the inner MATLAB function's logic and handles 
    the 'persistent fcount' via the module-level list `_compute_a_fcount`.
    """
    global _compute_a_fcount
    
    # Handle the persistent fcount initialization/reset call (if nargin == 0)
    if Ci is None:
        _compute_a_fcount[0] = 0
        return None, None
        
    fcount = _compute_a_fcount[0]

    if Type.lower() == 'c3':
        Vs = Vs_C3
        
        if np.any(g_m < np.inf):
            # With g_m (Mesophyll conductance): Vc and Ve solved from quadratic equations
            
            # Vc limited
            a_vc = 1.0 / g_m
            b_vc = -(MM_consts + Ci + (Rd + Vcmax) / g_m)
            c_vc = Vcmax * (Ci - Gamma_star + Rd / g_m)
            Vc = _sel_root(a_vc, b_vc, c_vc, -1) # -1 sign for the smallest root
            
            # Ve limited
            a_ve = 1.0 / g_m
            b_ve = -(Ci + 2 * Gamma_star + (Rd + Je * effcon) / g_m)
            c_ve = Je * effcon * (Ci - Gamma_star + Rd / g_m)
            Ve = _sel_root(a_ve, b_ve, c_ve, -1) # -1 sign for the smallest root
            
            CO2_per_electron = Ve / Je
        else:
            # Without g_m
            Vc = Vcmax * (Ci - Gamma_star) / (MM_consts + Ci)
            CO2_per_electron = (Ci - Gamma_star) / (Ci + 2 * Gamma_star) * effcon
            Ve = Je * CO2_per_electron
            
    else: # C4
        # C4 assumes no mesophyll resistance in this implementation
        Vc = Vcmax
        Vs = kpepcase * Ci
        CO2_per_electron = effcon 
        Ve = Je * CO2_per_electron

    # Smoothing min(Vc, Ve) -> V
    # sign(-Vc) logic is equivalent to sign(Gamma_star - Ci) in this context
    V = _sel_root(atheta, -(Vc + Ve), Vc * Ve, np.sign(Gamma_star - Ci)) 
    
    # Smoothing min(V, Vs) -> Ag (Gross Assimilation)
    Ag = _sel_root(0.98, -(V + Vs), V * Vs, -1)
    
    A = Ag - Rd
    fcount += 1
    _compute_a_fcount[0] = fcount
        
    biochem_out = {
        'A': A, 'Ag': Ag, 'Vc': Vc, 'Vs': Vs, 'Ve': Ve,
        'CO2_per_electron': CO2_per_electron, 'fcount': fcount
    }
    
    return A, biochem_out

def _temperature_function_c3(Tref, R, T, deltaHa):
    """Temperature function for C3 (exponential factor only)"""
    tempfunc1 = (1 - Tref / T)
    fTv = np.exp(deltaHa / (Tref * R) * tempfunc1)
    return fTv

def _high_temp_inhibtion_c3(Tref, R, T, deltaS, deltaHd):
    """High Temperature Inhibition Function for C3"""
    hightempfunc_num = (1 + np.exp((Tref * deltaS - deltaHd) / (Tref * R)))
    hightempfunc_deno = (1 + np.exp((deltaS * T - deltaHd) / (R * T)))
    fHTv = hightempfunc_num / hightempfunc_deno
    return fHTv

def _fluorescencemodel(ps, x, Kp, Kf, Kd, Knparams):
    """
    MATLAB's Fluorescencemodel translation.
    """
    Kno, alpha, beta = Knparams[0], Knparams[1], Knparams[2]
    
    # Kn = Kno * (1+beta).* x_alpha./(beta + x_alpha)
    x_alpha = np.power(x, alpha)
    Kn = Kno * (1 + beta) * x_alpha / (beta + x_alpha)
    
    fo0 = Kf / (Kf + Kp + Kd)        
    fo = Kf / (Kf + Kp + Kd + Kn)     
    fm = Kf / (Kf + Kd + Kn)          
    fm0 = Kf / (Kf + Kd)               
    fs = fm * (1 - ps)              
    
    # Handle division by zero
    eta = np.where(fo0 > 0, fs / fo0, 0.0)
    qQ = 1 - (fs - fo) / (fm - fo)
    qE = 1 - (fm - fo) / (fm0 - fo0)  
    
    return eta, qE, qQ, fs, fo, fm, fo0, fm0, Kn

def satvap(T_c):
    """Saturated vapour pressure (hPa) at temperature T (Celsius)"""
    return 6.107 * 10**(7.5 * T_c / (237.3 + T_c))

# --- Main Function ---

def biochemical(leafbio, meteo, options, constants, fV):
    """
    Translation of src/fluxes/biochemical.m
    """
    global _compute_a_fcount
    
    # --- Input Parsing and Setup ---
    
    tempcor = options['apply_T_corr']
    rhoa, Mair, R = constants['rhoa'], constants['Mair'], constants['R']
    Q = meteo['Q']
    Cs = meteo['Cs']
    T_celsius = meteo['T']
    T = meteo['T'] + 273.15 * (meteo['T'] < 200) # [K]
    eb, O, p = meteo['eb'], meteo['Oa'], meteo['p']
    
    Type = leafbio['Type']
    Vcmax25 = fV * leafbio['Vcmax25']
    BallBerrySlope, RdPerVcmax25, BallBerry0 = leafbio['BallBerrySlope'], leafbio['RdPerVcmax25'], leafbio['BallBerry0']
    Tref = 25 + 273.15
    
    Kc25, Ko25, spfy25 = 405, 279, 2444
    
    # Unit conversions (to bar)
    ppm2bar = 1e-6 * (p * 1e-3) 
    Cs = Cs * ppm2bar 
    O = (O * 1e-3) * (p * 1e-3) * (1 if Type.lower() == 'c3' else 0) 
    Kc25, Ko25 = Kc25 * 1e-6, Ko25 * 1e-3
    Gamma_star25 = 0.5 * O / spfy25      
    Rd25 = RdPerVcmax25 * Vcmax25
    
    effcon = 1/5 if Type.lower() == 'c3' else 1/6
    atheta = 0.8
    
    g_m = np.inf
    if 'g_m' in leafbio:
        g_m = leafbio['g_m'] * 1e6 
        
    stressfactor = leafbio['stressfactor']
    
    Knparams = np.array([leafbio['Kn0'], leafbio['Knalpha'], leafbio['Knbeta']])
    Kf, Kp = 0.05, 4.0
    Kd = np.maximum(0.8738, 0.0301 * T_celsius + 0.0773)

    # --- Temperature Corrections ---
    
    f = {'Vcmax': 1, 'Rd': 1, 'Kc': 1, 'Ko': 1, 'Gamma_star': 1}
    Ke = 1
    
    if tempcor:
        TDP = leafbio['TDP']
        if Type.lower() == 'c4':
            Q10, s1, s2, s3, s4, s5, s6 = TDP['Q10'], TDP['s1'], TDP['s2'], TDP['s3'], TDP['s4'], TDP['s5'], TDP['s6']
            
            fHTv_v = 1 + np.exp(s1 * (T - s2))
            fLTv_v = 1 + np.exp(s3 * (s4 - T))
            Vcmax = (Vcmax25 * np.power(Q10, 0.1 * (T - Tref))) / (fHTv_v * fLTv_v)
            
            fHTv_r = 1 + np.exp(s5 * (T - s6))
            Rd = (Rd25 * np.power(Q10, 0.1 * (T - Tref))) / fHTv_r
            
            Ke25 = 20000 * Vcmax25
            Ke = Ke25 * np.power(Q10, 0.1 * (T - Tref))
        
        else: # C3
            # Vcmax
            fTv_v = _temperature_function_c3(Tref, R, T, TDP['delHaV'])
            fHTv_v = _high_temp_inhibtion_c3(Tref, R, T, TDP['delSV'], TDP['delHdV'])
            f['Vcmax'] = fTv_v * fHTv_v
            
            # Rd
            fTv_r = _temperature_function_c3(Tref, R, T, TDP['delHaR'])
            fHTv_r = _high_temp_inhibtion_c3(Tref, R, T, TDP['delSR'], TDP['delHdR'])
            f['Rd'] = fTv_r * fHTv_r
            
            # Kc, Ko, Gamma_star (only TV factor)
            f['Kc'] = _temperature_function_c3(Tref, R, T, TDP['delHaKc'])
            f['Ko'] = _temperature_function_c3(Tref, R, T, TDP['delHaKo'])
            f['Gamma_star'] = _temperature_function_c3(Tref, R, T, TDP['delHaT'])
            
            Vcmax = Vcmax25 * f['Vcmax'] * stressfactor
            Rd = Rd25 * f['Rd'] * stressfactor
            Kc = Kc25 * f['Kc']
            Ko = Ko25 * f['Ko']
            
        Gamma_star = Gamma_star25 * f['Gamma_star']
        
    else: # No temperature correction
        Vcmax = Vcmax25 * stressfactor
        Rd = Rd25 * stressfactor
        Kc, Ko = Kc25, Ko25
        Gamma_star = Gamma_star25

    # --- Electron Transport (Je) ---
    po0 = Kp / (Kf + Kd + Kp)         
    Je = 0.5 * po0 * Q          

    # --- Photosynthesis constants ---
    if Type.lower() == 'c3':
        MM_consts = (Kc * (1 + O / Ko)) 
        Vs_C3 = (Vcmax / 2) 
        minCi = 0.3
    else:
        MM_consts = 0 
        Vs_C3 = 0     
        minCi = 0.1 

    # --- Ci Calculation (Iteration) ---
    
    RH = np.minimum(1, eb / satvap(T_celsius))

    # Reset persistent fcount before calculation
    _compute_a(None, None, None, None, None, None, None, None, None, None, None, None) 
    
    # Create the closure for computeA to be used by Ci_next
    computeA_fun = lambda x: _compute_a(x, Type, g_m, Vs_C3, MM_consts, Rd, Vcmax, Gamma_star, Je, effcon, atheta, Ke)

    if np.all(BallBerry0 == 0):
        # Non-iterative: A = 0 assumption
        Ci, _ = _ball_berry(Cs, RH, None, BallBerrySlope, BallBerry0, minCi)
    else:
        # Iterative solution: Ci = f(A(Ci))
        Ci_next_fun = lambda x: _ci_next(x, Cs, RH, minCi, BallBerrySlope, BallBerry0, computeA_fun, ppm2bar)
        tol = 1e-7
        # fixedp_brent_ari requires a function that returns (error, Ci_next_val)
        Ci = _fixedp_brent_ari(Ci_next_fun, Cs, None, tol) 

    # Final calculation of A
    A, biochem_out = computeA_fun(Ci)

    # --- Fluxes and Resistances ---
    Ag = biochem_out['Ag']
    CO2_per_electron = biochem_out['CO2_per_electron']

    gs = np.maximum(0, 1.6 * A * ppm2bar / (Cs - Ci))     
    Ja = Ag / CO2_per_electron   
    rcw = (rhoa / (Mair * 1e-3)) / gs
    
    # --- Fluorescence ---
    
    ps = po0 * Ja / Je               
    nan_ps = np.isnan(ps)
    if np.any(nan_ps):
        ps = np.where(nan_ps, np.broadcast_to(po0, ps.shape)[nan_ps], ps)
    
    ps_rel = np.maximum(0, 1 - ps / po0)
    
    eta, qE, qQ, fs, fo, fm, fo0, fm0, Kn = _fluorescencemodel(ps, ps_rel, Kp, Kf, Kd, Knparams)
    Kpa = ps / fs * Kf
    
    # --- Convert back to ppm ---
    
    Cc = None
    if g_m < np.inf:
        Cc_bar = Ci - A / g_m
        Cc = Cc_bar / ppm2bar
        
    Ci_ppm = Ci / ppm2bar
    
    # --- Collect outputs ---
    
    biochem_out['A'] = A
    biochem_out['Ci'] = Ci_ppm
    if Cc is not None:
        biochem_out['Cc'] = Cc
    biochem_out['rcw'] = rcw
    biochem_out['gs'] = gs
    biochem_out['RH'] = RH
    biochem_out['Vcmax'] = Vcmax
    biochem_out['Rd'] = Rd
    biochem_out['Ja'] = Ja
    biochem_out['ps'] = ps
    biochem_out['ps_rel'] = ps_rel
    biochem_out['Kd'] = Kd
    biochem_out['Kn'] = Kn
    biochem_out['NPQ'] = Kn / (Kf + Kd) 
    biochem_out['Kf'] = Kf
    biochem_out['Kp0'] = Kp
    biochem_out['Kp'] = Kpa
    biochem_out['eta'] = eta
    biochem_out['qE'] = qE
    biochem_out['fs'] = fs
    biochem_out['ft'] = fs
    biochem_out['SIF'] = fs * Q
    biochem_out['fo0'] = fo0
    biochem_out['fm0'] = fm0
    biochem_out['fo'] = fo
    biochem_out['fm'] = fm
    biochem_out['Fm_Fo'] = fm / fo
    biochem_out['Ft_Fo'] = fs / fo
    biochem_out['qQ'] = qQ
    biochem_out['Phi_N'] = Kn / (Kn + Kp + Kf + Kd)
    
    return biochem_out