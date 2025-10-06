import numpy as np
from scipy.special import exp1 # Exponential Integral E1

# --- Helper Functions (Private to Module) ---

def _calctav(alfa, nr):
    """Translation of MATLAB's calctav function."""
    rd = np.pi / 180
    n2 = nr**2
    np_val = n2 + 1
    nm = n2 - 1
    a_val = (nr + 1)**2 / 2
    k_val = -(n2 - 1)**2 / 4
    sa = np.sin(alfa * rd)

    # Use where to handle the scalar alfa!=90 condition
    b2 = sa**2 - np_val / 2
    
    b1_mask = (alfa != 90)
    b1 = np.where(b1_mask, np.sqrt((b2**2 + k_val)), 0)
    
    b = b1 - b2
    b3, a3 = b**3, a_val**3
    
    # Suppress warnings for log/div by zero which are handled by masks
    with np.errstate(divide='ignore', invalid='ignore'):
        
        ts = (k_val**2 / (6 * b3) + k_val / b - b / 2) - \
             (k_val**2 / (6 * a3) + k_val / a_val - a_val / 2)

        tp1 = -2 * n2 * (b - a_val) / (np_val**2)
        
        # Mask for nm = 0 (n=1)
        nm_mask = (nm == 0)
        
        tp2_safe = np.where(nm_mask, 0, -2 * n2 * np_val * np.log(b / a_val) / (nm**2))
        tp3_safe = np.where(nm_mask, 0, n2 * (1 / b - 1 / a_val) / 2)
        
        tp4_num = (2 * np_val * b - nm**2) / (2 * np_val * a_val - nm**2)
        tp4_safe = np.where(nm_mask, 0, 16 * n2**2 * (n2**2 + 1) * np.log(tp4_num) / (np_val**3 * nm**2))
        
        tp5_den_b = 2 * np_val * b - nm**2
        tp5_den_a = 2 * np_val * a_val - nm**2
        tp5_safe = np.where(nm_mask, 0, 16 * n2**3 * (1 / tp5_den_b - 1 / tp5_den_a) / (np_val**3))
        
        tp = tp1 + tp2_safe + tp3_safe + tp4_safe + tp5_safe
        
        # tav is only computed if sa != 0
        sa_mask = (sa != 0)
        tav = np.where(sa_mask, (ts + tp) / (2 * sa**2), 4 * nr / ((nr + 1)**2))
        
    return tav

# --- Main Function ---

def fluspect_B_CX(spectral, leafbio, optipar):
    """
    Translation of src/RTMs/fluspect_B_CX.m
    """
    
    # --- parameters ---
    ndub = 15
    int_factor = 5 # Used as spectral spacing reduction factor for fluorescence
    
    Cab, Cca, V2Z, Cw, Cdm, Cs = leafbio['Cab'], leafbio['Cca'], leafbio['V2Z'], leafbio['Cw'], leafbio['Cdm'], leafbio['Cs']
    Cant, Cbc, Cp, N, fqe = leafbio['Cant'], leafbio['Cbc'], leafbio['Cp'], leafbio['N'], leafbio['fqe']

    nr, Kdm, Kab = optipar['nr'], optipar['Kdm'], optipar['Kab']
    
    if V2Z == -999:
        Kca = optipar['Kca']
    else:
        Kca = (1 - V2Z) * optipar['KcaV'] + V2Z * optipar['KcaZ']
        
    Kw, Ks, Kant = optipar['Kw'], optipar['Ks'], optipar['Kant']
    Kp = optipar.get('Kp', np.zeros_like(Kab))
    Kcbc = optipar.get('Kcbc', np.zeros_like(Kab))
    phi = optipar['phi']

    # --- PROSPECT calculations ---
    # Kall = (Cab*Kab + Cca*Kca + Cdm*Kdm + Cw*Kw + Cs*Ks + Cant*Kant + Cp*Kp + Cbc*Kcbc)/N;
    Kall = (Cab * Kab + Cca * Kca + Cdm * Kdm + Cw * Kw + Cs * Ks + Cant * Kant + Cp * Kp + Cbc * Kcbc) / N
    
    # Non-conservative scattering
    j = (Kall > 0)
    
    # Taylor series approximation of integral (t1 + t2) for small Kall
    t1 = (1 - Kall) * np.exp(-Kall)
    t2 = Kall**2 * exp1(Kall) # exp1 is the scipy.special.exp1 function (exponential integral)
    
    tau = np.ones_like(t1)
    tau[j] = t1[j] + t2[j]
    
    # Specific absorption coefficients (used for Rnsun_Cab in RTMo)
    kChlrel, kCarrel = np.zeros_like(t1), np.zeros_like(t1)
    kChlrel[j] = Cab * Kab[j] / (Kall[j] * N)
    kCarrel[j] = Cca * Kca[j] / (Kall[j] * N)

    # Interface reflection/transmission (59 deg is standard for PROSPECT)
    talf = _calctav(59, nr)
    ralf = 1 - talf
    t12 = _calctav(90, nr) # T(90, n_water -> n_air)
    r12 = 1 - t12
    t21 = t12 / (nr**2) # T(90, n_air -> n_water)
    r21 = 1 - t21

    # Top surface side (top interface + mesophyll stack reflection/transmission)
    denom = 1 - r21 * r21 * tau**2
    Ta = talf * tau * t21 / denom
    Ra = ralf + r21 * tau * Ta

    # Bottom surface side (mesophyll stack reflection/transmission only)
    t = t12 * tau * t21 / denom
    r = r12 + r21 * t
    
    # --- Stokes equations for next N-1 layers ---
    
    D = np.sqrt((1 + r + t) * (1 + r - t) * (1 - r + t) * (1 - r - t))
    rq, tq = r**2, t**2
    a = (1 + rq - tq + D) / (2 * r)
    b = (1 - rq + tq + D) / (2 * t)

    # Special case: N=1 (a=Inf, b=1) or a+r=1 (r+t >= 1, zero absorption)
    # General case for N-1 layers
    bNm1 = b**(N - 1)
    bN2 = bNm1**2
    a2 = a**2
    denom_stokes = a2 * bN2 - 1
    
    Rsub = a * (bN2 - 1) / denom_stokes
    Tsub = bNm1 * (a2 - 1) / denom_stokes
    
    # Case of zero absorption (r+t >= 1)
    j_zero_abs = (r + t >= 1)
    if np.any(j_zero_abs):
        Tsub[j_zero_abs] = t[j_zero_abs] / (t[j_zero_abs] + (1 - t[j_zero_abs]) * (N - 1))
        Rsub[j_zero_abs] = 1 - Tsub[j_zero_abs]
        
    # Combine top layer (Ra, Ta) with N-1 stack (Rsub, Tsub)
    denom_final = 1 - Rsub * r
    tran = Ta * Tsub / denom_final
    refl = Ra + Ta * Rsub * t / denom_final

    leafopt = {'refl': refl, 'tran': tran, 'kChlrel': kChlrel, 'kCarrel': kCarrel}

    # --- Fluorescence Doubling Method Setup ---

    # Remove the interfaces to get mesophyll-only properties (rho, tau)
    Rb = (refl - ralf) / (talf * t21 + (refl - ralf) * r21)
    Z = tran * (1 - Rb * r21) / (talf * t21)
    
    rho_meso = (Rb - r21 * Z**2) / (1 - (r21 * Z)**2)
    tau_meso = (1 - Rb * r21) / (1 - (r21 * Z)**2) * Z
    
    rho_meso = np.maximum(rho_meso, 0)
    
    # Derive Kubelka-Munk s and k (for single elementary layer, not the N-layer stack)
    r, t = rho_meso, tau_meso
    I_rt = (r + t) < 1
    
    # D, a, b re-calculated for mesophyll layer properties
    D = np.ones_like(r) * np.inf # Placeholder for non-conservative scattering
    D[I_rt] = np.sqrt((1 + r[I_rt] + t[I_rt]) * (1 + r[I_rt] - t[I_rt]) * (1 - r[I_rt] + t[I_rt]) * (1 - r[I_rt] - t[I_rt]))
    a = np.ones_like(r) * np.inf
    b = np.ones_like(r) * np.inf
    a[I_rt] = (1 + r[I_rt]**2 - t[I_rt]**2 + D[I_rt]) / (2 * r[I_rt])
    b[I_rt] = (1 - r[I_rt]**2 + t[I_rt]**2 + D[I_rt]) / (2 * t[I_rt])

    s = r / t # Default for r+t=1 (conservative scattering)
    k = np.log(b)
    
    # Non-conservative scattering (a>1)
    I_a = (a > 1) & (a != np.inf)
    s[I_a] = 2 * a[I_a] / (a[I_a]**2 - 1) * np.log(b[I_a])
    k[I_a] = (a[I_a] - 1) / (a[I_a] + 1) * np.log(b[I_a])
    
    kChl = kChlrel * k

    # --- Fluorescence of the leaf mesophyll layer ---

    if fqe > 0:
        # Interpolate K, s, kChl, r21, tau, etc., to the excitation/fluorescence grids
        wle = np.arange(400, 751, int_factor)
        wlf = np.arange(640, 851, 4)
        
        # Use np.interp (faster, linear) or _interp1 (scipy, spline)
        interp_func = np.interp
        
        # Excitation wavelengths
        k_iwle = interp_func(wle, spectral['wlP'], k)
        s_iwle = interp_func(wle, spectral['wlP'], s)
        kChl_iwle = interp_func(wle, spectral['wlP'], kChl)
        r21_iwle = interp_func(wle, spectral['wlP'], r21)
        rho_iwle = interp_func(wle, spectral['wlP'], rho_meso)
        tau_iwle = interp_func(wle, spectral['wlP'], tau_meso)
        talf_iwle = interp_func(wle, spectral['wlP'], talf)
        
        # Fluorescence wavelengths
        Iwlf = np.isin(spectral['wlP'], wlf)
        
        k_iwlf = interp_func(wlf, spectral['wlP'], k)
        s_iwlf = interp_func(wlf, spectral['wlP'], s)
        r21_iwlf = interp_func(wlf, spectral['wlP'], r21)
        rho_iwlf = interp_func(wlf, spectral['wlP'], rho_meso)
        tau_iwlf = interp_func(wlf, spectral['wlP'], tau_meso)
        
        eps = 2**(-ndub)

        # Initial conditions for doubling
        te = 1 - (k_iwle + s_iwle) * eps
        tf = 1 - (k_iwlf + s_iwlf) * eps
        re = s_iwle * eps
        rf = s_iwlf * eps
        
        # Sigmoid matrix [n_wlf, n_wle]
        wlf_exp = wlf[:, np.newaxis]
        wle_exp = wle[np.newaxis, :]
        sigmoid = 1 / (1 + np.exp(-wlf_exp / 10) * np.exp(wle_exp / 10))
        
        # Initial fluorescence matrices: [n_wlf, n_wle]
        Mb = int_factor * fqe * (0.5 * phi[Iwlf])[:, np.newaxis] * eps * kChl_iwle[np.newaxis, :] * sigmoid
        Mf = Mb.copy()
        
        # Doubling routine
        for i in range(ndub):
            # Element-wise operations
            xe = te / (1 - re * re); ten = te * xe; ren = re * (1 + ten);
            xf = tf / (1 - rf * rf); tfn = tf * xf; rfn = rf * (1 + tfn);
            
            # Matrix-vector outer products (A11, A12, A21, A22 are expanded)
            xf_col, xe_row = xf[:, np.newaxis], xe[np.newaxis, :]
            rf_col, re_row = rf[:, np.newaxis], re[np.newaxis, :]
            
            # A terms for Mf/Mb update
            A11 = xf_col + xe_row # [n_wlf, n_wle]
            A12 = (xf_col * xe_row) * (rf_col + re_row)
            A21 = 1 + (xf_col * xe_row) * (1 + rf_col * re_row)
            A22 = (xf_col * rf_col) + (xe_row * re_row)
            
            Mfn = Mf * A11 + Mb * A12
            Mbn = Mb * A21 + Mf * A22
            
            te, re, tf, rf = ten, ren, tfn, rfn
            Mf, Mb = Mfn, Mbn
            
        # --- Add interfaces back ---
        
        # Rb, Z for the whole leaf (using mesophyll properties to link to interfaces)
        Rb_whole = refl
        
        # Interpolate terms back to fluorescence grid
        Rb_iwle = interp_func(wle, spectral['wlP'], Rb_whole)
        
        Xe = talf_iwle[np.newaxis, :] / (1 - r21_iwle * Rb_iwle)[np.newaxis, :] # [1, n_wle]
        Xf = t21[Iwlf][:, np.newaxis] / (1 - r21[Iwlf] * Rb_whole[Iwlf])[:, np.newaxis] # [n_wlf, 1]
        
        Ye = tau_iwle[np.newaxis, :] * r21_iwle[np.newaxis, :] / (1 - rho_iwle * r21_iwle)[np.newaxis, :] # [1, n_wle]
        Yf = tau_iwlf[:, np.newaxis] * r21_iwlf[:, np.newaxis] / (1 - rho_iwlf * r21_iwlf)[:, np.newaxis] # [n_wlf, 1]

        # Final fluorescence matrices (gn, fn are Mb, Mf)
        A = Xe * (1 + Ye * Yf) * Xf
        B = Xe * (Ye + Yf) * Xf
        
        gn = A * Mb + B * Mf
        fn = A * Mf + B * Mb
        
        leafopt['Mb'] = gn
        leafopt['Mf'] = fn
        
    return leafopt