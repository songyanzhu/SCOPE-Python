import numpy as np

# --- Helper Functions (As defined in RTMt_planck.py) ---
def _stefan_boltzmann(T_C, constants):
    sigmaSB = constants['sigmaSB']
    C2K = constants['C2K']
    H = sigmaSB * (T_C + C2K)**4
    return H

# --- Main Function ---

def RTMt_sb(constants, rad, soil, leafbio, canopy, gap, Tcu, Tch, Tsu, Tsh, obsdir=False, spectral=None):
    """
    Translation of src/RTMs/RTMt_sb.m
    Simplified Thermal RTM using Stefan-Boltzmann (non-spectral).
    """
    
    # --- 0.1 parameters ---
    
    nl = canopy['nlayers']
    lidf = canopy['lidf']
    Ps = gap['Ps']
    
    # Thermal properties (assumed to be scalar from leafbio/soil)
    rho = leafbio['rho_thermal']
    tau = leafbio['tau_thermal']
    rs = soil['rs_thermal']
    epsc = 1 - rho - tau
    epss = 1 - rs
    
    LAI = canopy['LAI']
    dx = 1 / nl
    iLAI = LAI * dx
    
    # RTM coefficients (assumed to be a scalar/vector slice for the thermal band)
    # The MATLAB code takes the *last* element, often used for thermal band (e.g., 10um)
    
    # Take the last wavelength's coefficients as representative of the thermal band
    last_idx = -1
    Xdd = rad['Xdd'][:, last_idx]
    R_dd = rad['R_dd'][:, last_idx]
    rho_dd = rad['rho_dd'][:, last_idx]
    tau_dd = rad['tau_dd'][:, last_idx]
    
    Xss = np.tile(rad['Xss'][:, last_idx], nl)
    R_sd = rad['R_sd'][:, last_idx]
    Xsd = rad['Xsd'][:, last_idx]
    
    # --- 1. calculation of upward and downward fluxes ---

    # 1.1 Radiance by components (H = epsilon * sigma * T^4)
    Hcsu3 = epsc * _stefan_boltzmann(Tcu, constants)
    Hcsh = epsc * _stefan_boltzmann(Tch, constants)
    Hssu = epss * _stefan_boltzmann(Tsu, constants)
    Hssh = epss * _stefan_boltzmann(Tsh, constants)

    # 1.2 Averaging Hcsu3 (Hemispherical emittance by leaf layers)
    if Hcsu3.ndim > 1:
        # Average over angles (13x36 matrix) for each layer
        Hcsu_angled = Hcsu3.reshape(nlinc, nlazi, nl)
        # Weight angles by lidf (approximation: dot product then reshape)
        Hcsu_mean = np.sum(Hcsu_angled * lidf[:, np.newaxis, np.newaxis] / np.sum(lidf), axis=(0, 1))
        Hcsu = Hcsu_mean
    else:
        Hcsu = Hcsu3
    
    Hc = Hcsu * Ps[:nl] + Hcsh * (1 - Ps[:nl]) # [nl]
    Hs = Hssu * Ps[nl] + Hssh * (1 - Ps[nl])   # [1]

    # 1.3 Diffuse radiation (Scalar RTM, vectorized over nl for leaf layers)
    U, Emin, Eplu = np.zeros(nl + 1), np.zeros(nl + 1), np.zeros(nl + 1)
    
    U[nl] = Hs
    
    # Downward RTM (Bottom to Top)
    Y = np.zeros(nl)
    for j in range(nl - 1, -1, -1):
        Y[j] = (rho_dd[j] * U[j+1] + Hc[j] * iLAI) / (1 - rho_dd[j] * R_dd[j+1])
        U[j] = tau_dd[j] * (R_dd[j+1] * Y[j] + U[j+1]) + Hc[j] * iLAI

    # Upward RTM (Top to Bottom)
    for j in range(nl):
        # Es_(j) is 0 for thermal RTM
        Emin[j+1] = Xsd[j] * 0 + Xdd[j] * Emin[j] + Y[j]
        Eplu[j] = R_sd[j] * 0 + R_dd[j] * Emin[j] + U[j]
        
    Eplu[nl] = R_sd[nl - 1] * 0 + R_dd[nl - 1] * Emin[nl] + Hs
    Eoutte = Eplu[0] # TOC Hemispherical Emitted Flux

    # 1.4 Directional radiation (if obsdir is True)
    piLot = None
    if obsdir:
        K = gap['K']
        vb, vf = rad['vb'][last_idx], rad['vf'][last_idx]
        
        # Directional Emitted/Scattered from Vegetation (piLov) and Soil (piLos)
        piLov = iLAI * (K * Hcsh @ (gap['Po'][:nl] - gap['Pso'][:nl]) + 
                        K * Hcsu @ gap['Pso'][:nl] +
                        (vb * Emin[:nl] + vf * Eplu[1:]).T @ gap['Po'][:nl])
        
        piLos = (Hssh * (gap['Po'][nl] - gap['Pso'][nl]) + Hssu * gap['Pso'][nl])
        piLot = piLov + piLos
        
    # --- 2. total net fluxes ---
    
    # Rnuc = epsc*(Emin(1:end-1) + Eplu(2:end)) - 2*(Hcsu)
    Rnuc = epsc * (Emin[:nl] + Eplu[1:]) - 2 * Hcsu # [nl]
    Rnhc = epsc * (Emin[:nl] + Eplu[1:]) - 2 * Hcsh # [nl]
    Rnus = epss * (Emin[nl] - Hssu) # [1]
    Rnhs = epss * (Emin[nl] - Hssh) # [1]

    # --- 3. Write the output to the rad structure ---
    rad['Emint'] = Emin
    rad['Eplut'] = Eplu
    rad['Eoutte'] = Eoutte
    
    # Net thermal fluxes (used in energy balance iteration)
    rad['Rnuct'] = Rnuc
    rad['Rnhct'] = Rnhc
    rad['Rnust'] = Rnus
    rad['Rnhst'] = Rnhs
    
    if obsdir:
        rad['Lote'] = piLot / np.pi
    
    return rad