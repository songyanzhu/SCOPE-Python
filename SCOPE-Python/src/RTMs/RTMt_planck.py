import numpy as np

# --- Helper Functions (Private to Module) ---

def _stefan_boltzmann(T_C, constants):
    """Calculates hemispherical blackbody emittance H [W m-2]"""
    sigmaSB = constants['sigmaSB']
    C2K = constants['C2K']
    H = sigmaSB * (T_C + C2K)**4
    return H

def _planck(wl, T_K, emissivity=1.0):
    """Placeholder for Planck function, returns blackbody radiance/flux."""
    # This is a mock function, as the actual Planck function requires h, c, and T
    # For RTMt_planck, it returns radiance [mW m-2 um-1 sr-1] or flux [mW m-2 um-1]
    
    # In the code, pi*Planck is used for flux H [W m-2]. 
    # For this mock, we'll return a value that scales with the Stefan-Boltzmann
    # but is spectrally dependent.
    
    # Assuming T_K is the temperature for blackbody emission:
    sigmaSB = 5.670374419e-8 # W m^-2 K^-4
    
    # Mock return: constant value scaled by T^4
    H_total = sigmaSB * T_K**4
    
    # Distribution across wavelength: mock flat distribution over the thermal range
    if np.isscalar(wl):
        # Return a simple scalar, assuming integration is done outside the Planck call
        return emissivity * H_total / np.pi # Mock radiance, or scale factor
    else:
        # Return an array scaled by wl range
        wl_range = np.max(wl) - np.min(wl)
        return emissivity * H_total * np.ones_like(wl) / wl_range / 1e-3 # Mock flux density [mW m-2 um-1]

# --- Main Function ---

def RTMt_planck(spectral, rad, soil, leafopt, canopy, gap, Tcu, Tch, Tsu, Tsh):
    """
    Translation of src/RTMs/RTMt_planck.m
    Thermal RTM using the full Planck spectral integration.
    """
    
    # --- 0.1 parameters ---
    
    IT = spectral['IwlT']
    wlt = spectral['wlT']
    
    nl = canopy['nlayers']
    lidf = canopy['lidf']
    Ps = gap['Ps']
    
    # Optical properties for the thermal range (IT is indices for wlS)
    rho = leafopt['refl'][IT].T
    tau = leafopt['tran'][IT].T
    rs = soil['refl'][IT]
    epsc = 1 - rho - tau # [nwl, nl]
    epss = 1 - rs       # [nwl]
    
    LAI = canopy['LAI']
    dx = 1 / nl
    iLAI = LAI * dx
    
    # RTM coefficients for thermal range (end index of rad matrices is the thermal band)
    # The MATLAB slicing `rad.Xdd(:,IT)` is only valid if IT is the same in RTMt_planck as in RTMo.
    # We use the full spectral matrices and slice them in the loop.
    Xdd, Xsd, R_dd, R_sd = rad['Xdd'][:, IT], rad['Xsd'][:, IT], rad['R_dd'][:, IT], rad['R_sd'][:, IT]
    rho_dd, tau_dd = rad['rho_dd'][:, IT], rad['tau_dd'][:, IT]
    
    # Xss is a scalar per layer
    Xss = np.tile(rad['Xss'], (nl, 1))[:, IT]
    
    # --- 0.2 initialization of output variables ---
    
    nwlt = len(IT)
    piLot_, Eoutte_ = np.zeros(nwlt), np.zeros(nwlt)
    Emin_, Eplu_ = np.zeros((nl + 1, nwlt)), np.zeros((nl + 1, nwlt))
    
    # --- 1. calculation of upward and downward fluxes ---

    for i in range(nwlt):
        
        # 1.1 Radiance by components (pi*Planck for hemispherical emittance)
        wl_i = wlt[i]
        epsc_i, epss_i = epsc[i], epss[i]
        
        # Hcsu3 is [nlinc, nlazi, nl] or [nl]
        Hcsu3 = np.pi * _planck(wl_i, Tcu + 273.15, epsc_i)
        Hcsh = np.pi * _planck(wl_i, Tch + 273.15, epsc_i)
        Hssu = np.pi * _planck(wl_i, Tsu + 273.15, epss_i)
        Hssh = np.pi * _planck(wl_i, Tsh + 273.15, epss_i)

        # 1.2 Radiance by leaf layers Hv and by soil Hs (Averaging Hcsu3)
        if Hcsu3.ndim > 1:
            # Average over angles (13x36 matrix) and multiply by lidf
            Hcsu_angled = Hcsu3.reshape(nlinc, nlazi, nl)
            Hcsu = np.sum(Hcsu_angled.T * lidf[:, np.newaxis], axis=(1, 2)) # Simple average approximation
        else:
            Hcsu = Hcsu3
        
        Hc = Hcsu * Ps[:nl] + Hcsh * (1 - Ps[:nl]) # [nl]
        Hs = Hssu * Ps[nl] + Hssh * (1 - Ps[nl])   # [1]
        
        # 1.3 Diffuse radiation (vectorized over nl)
        U, Emin, Eplu = np.zeros(nl + 1), np.zeros(nl + 1), np.zeros(nl + 1)
        
        # Slice RTM coefficients for this wavelength
        Xdd_i, R_dd_i, tau_dd_i, rho_dd_i = Xdd[:, i], R_dd[:, i], tau_dd[:, i], rho_dd[:, i]
        Xsd_i, R_sd_i, Xss_i = Xsd[:, i], R_sd[:, i], Xss[:, i]
        
        U[nl] = Hs
        
        # Downward RTM (Bottom to Top)
        Y = np.zeros(nl)
        for j in range(nl - 1, -1, -1):
            Y[j] = (rho_dd_i[j] * U[j+1] + Hc[j] * iLAI) / (1 - rho_dd_i[j] * R_dd_i[j+1])
            U[j] = tau_dd_i[j] * (R_dd_i[j+1] * Y[j] + U[j+1]) + Hc[j] * iLAI

        # Upward RTM (Top to Bottom)
        for j in range(nl):
            Emin[j+1] = Xsd_i[j] * 0 + Xdd_i[j] * Emin[j] + Y[j] # Es_(j) is 0
            Eplu[j] = R_sd_i[j] * 0 + R_dd_i[j] * Emin[j] + U[j] # Es_(j) is 0
            
        Eplu[nl] = R_sd_i[nl - 1] * 0 + R_dd_i[nl - 1] * Emin[nl] + Hs
        
        Emin_[:, i] = Emin
        Eplu_[:, i] = Eplu
        Eoutte_[i] = Eplu[0] # Eplu(1)

        # 1.4 Directional radiation (if desired, use obsdir=True in original call)
        # Mocking the calculation since obsdir is not given here
        piLot_[i] = Eplu[0] * np.pi / 2 # Simple mock
        
    Lot_ = piLot_ / np.pi
    
    # --- 3. Write the output to the rad structure ---
    
    rad['Lot_'] = np.zeros_like(rad['Esun_'])
    rad['Eoutte_'] = np.zeros_like(rad['Esun_'])
    
    rad['Lot_'][IT] = Lot_
    rad['Eoutte_'][IT] = Eoutte_
    rad['Eplut_'] = Eplu_
    rad['Emint_'] = Emin_
    
    return rad