import numpy as np

# --- Helper Functions (Private to Module) ---

def _psim(z, L, unst, st, x):
    """Stability correction function for momentum (Paulson, 1970). z is z-d."""
    pm = np.zeros_like(z)
    
    if np.any(unst):
        # Unstable: pm = 2*log((1+x)/2)+log((1+x.^2)/2) -2*atan(x)+pi/2
        pm = np.where(unst, 
                      (2 * np.log((1 + x) / 2) + np.log((1 + x**2) / 2) - 2 * np.arctan(x) + np.pi / 2),
                      pm)
    
    if np.any(st):
        # Stable: pm = -5*z./L
        pm = np.where(st, -5 * z / L, pm)
        
    return pm

def _psih(z, L, unst, st, x):
    """Stability correction function for heat (Paulson, 1970). z is z-d."""
    ph = np.zeros_like(z)
    
    if np.any(unst):
        # Unstable: ph = 2*log((1+x.^2)/2)
        ph = np.where(unst, 2 * np.log((1 + x**2) / 2), ph)
        
    if np.any(st):
        # Stable: ph = -5*z./L
        ph = np.where(st, -5 * z / L, ph)
        
    return ph

def _phstar(z, zR, d, L, st, unst, x):
    """Stability correction function for the canopy layer (Paulson, 1970)."""
    phs = np.zeros_like(z)
    
    if np.any(unst):
        # Unstable: phs = (z-d)/(zR-d)*(x.^2-1)./(x.^2+1)
        phs = np.where(unst, (z - d) / (zR - d) * (x**2 - 1) / (x**2 + 1), phs)
        
    if np.any(st):
        # Stable: phs = -5*z./L
        phs = np.where(st, -5 * z / L, phs)
        
    return phs

# --- Main Function ---

def resistances(constants, soil, canopy, meteo):
    """
    Translation of src/fluxes/resistances.m
    """
    
    # --- Parameters ---
    kappa, Cd, LAI, rwc = constants['kappa'], canopy['Cd'], canopy['LAI'], canopy['rwc']
    z0m, d, h = canopy['zo'], canopy['d'], canopy['hc']
    z, u, L = meteo['z'], np.maximum(0.3, meteo['u']), meteo['L']
    rbs = soil['rbs']
    
    zr = 2.5 * h 
    n = Cd * LAI / (2 * kappa**2) 

    # --- Stability Correction ---
    unst = (L < 0) & (L > -500)
    st = (L > 0) & (L < 500)
    
    x = np.where(unst, (1 - 16 * z / L)**(1 / 4), 1) 

    pm_z = _psim(z - d, L, unst, st, x)
    pm_h = _psim(h - d, L, unst, st, x)
    ph_z = _psih(z - d, L, unst, st, x)
    
    ph_zr_at_zr = _psih(zr - d, L, unst, st, x)
    ph_zr = np.where(z >= zr, ph_zr_at_zr, ph_z)
    
    phs_zr = _phstar(zr, zr, d, L, st, unst, x)
    phs_h = _phstar(h, zr, d, L, st, unst, x)

    # Friction velocity (ustar)
    ustar = np.maximum(0.001, kappa * u / (np.log((z - d) / z0m) - pm_z)) 
    
    # Kinematic diffusivity (Kh)
    Kh_zr = kappa * ustar * (zr - d)                  
    Kh = Kh_zr
    Kh = np.where(unst, Kh_zr * (1 - 16 * (h - d) / L)**0.5, Kh)
    Kh = np.where(st, Kh_zr * (1 + 5 * (h - d) / L)**-1, Kh)

    # Wind speed at h (uh) and z0m (uz0)
    uh = np.maximum(ustar / kappa * (np.log((h - d) / z0m) - pm_h), 0.01)
    uz0 = uh * np.exp(n * ((z0m + d) / h - 1))                      

    # --- Resistances ---

    # rai (Inertial sublayer)
    rai = np.where(z > zr, 
                   1 / (kappa * ustar) * (np.log((z - d) / (zr - d)) - ph_z + ph_zr),
                   0)
    
    # rar (Roughness sublayer)
    rar = 1 / (kappa * ustar) * ((zr - h) / (zr - d) - phs_zr + phs_h)

    # rac (Canopy layer)
    log_term = (np.log((np.exp(n) - 1) / (np.exp(n) + 1)) - 
                np.log((np.exp(n * (z0m + d) / h) - 1) / (np.exp(n * (z0m + d) / h) + 1)))
    rac = h * np.sinh(n) / (n * Kh) * log_term
    
    # rws (Within canopy, soil)
    log_term_ws = (np.log((np.exp(n * (z0m + d) / h) - 1) / (np.exp(n * (z0m + d) / h) + 1)) - 
                   np.log((np.exp(n * (0.01) / h) - 1) / (np.exp(n * (0.01) / h) + 1)))
    rws = h * np.sinh(n) / (n * Kh) * log_term_ws

    raa = rai + rar + rac
    rawc = rwc
    raws = rws + rbs

    # --- Collect Outputs ---
    resist_out = {
        'ustar': ustar, 'raa': raa, 'rawc': rawc, 'raws': raws,
        'rai': rai, 'rar': rar, 'rac': rac, 'rws': rws, 'Kh': Kh, 'uz0': uz0
    }

    return resist_out