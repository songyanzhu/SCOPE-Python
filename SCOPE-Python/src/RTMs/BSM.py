import numpy as np

# --- Helper Function ---

def _tav(alfa, nr):
    """
    Translation of MATLAB's tav function, used in soilwat.
    Calculates the average escape probability for a beam of light.
    """
    n2 = nr**2
    np_val = n2 + 1
    nm = n2 - 1
    
    deg2rad = np.pi / 180
    sin_a = np.sin(alfa * deg2rad)
    
    if alfa != 0:
        B2 = sin_a**2 - np_val / 2
        
        # Calculate k (constant related to refractive index)
        k = -(n2 - 1)**2 / 4
        
        # B1 calculation handles alfa=90 case implicitly via the (alfa!=90) term in MATLAB, 
        # but in Python with numpy, we use where to handle the division by zero if present, 
        # though B1 itself is part of a square root. The MATLAB implementation is slightly 
        # simplified by using a boolean product: (alfa~=90) * sqrt(...)
        # We need to handle the potential complexity for alpha=90 where B1=0
        if alfa == 90:
            B1 = np.zeros_like(nr)
        else:
            B1 = np.sqrt(B2**2 + k)
        
        b = B1 - B2
        a = (nr + 1)**2 / 2
        
        b3 = b**3
        a3 = a**3
        
        # Suppress invalid value warnings for log/div by zero which are handled by where/masks
        with np.errstate(divide='ignore', invalid='ignore'):
            
            # ts (Transmittance S-polarization component)
            ts = (k**2 / (6 * b3) + k / b - b / 2) - \
                 (k**2 / (6 * a3) + k / a - a / 2)
                 
            # tp (Transmittance P-polarization component)
            tp1 = -2 * n2 * (b - a) / (np_val**2)
            tp2 = -2 * n2 * np_val * np.log(b / a) / (nm**2)
            tp3 = n2 * (1 / b - 1 / a) / 2
            
            # tp4 and tp5 need masking for nm=0 (n=1)
            nm_mask = (nm == 0)
            
            tp4 = np.zeros_like(nr)
            tp5 = np.zeros_like(nr)
            
            # Apply calculation only where nm != 0
            safe_nm = np.where(nm_mask, 1, nm)
            
            tp4_calc = 16 * n2**2 * (n2**2 + 1) * np.log((2 * np_val * b - nm**2) / (2 * np_val * a - nm**2)) / (np_val**3 * safe_nm**2)
            tp5_calc = 16 * n2**2 * n2 * (1 / (2 * np_val * b - nm**2) - 1 / (2 * np_val * a - nm**2)) / (np_val**3)
            
            tp4 = np.where(nm_mask, 0, tp4_calc)
            tp5 = np.where(nm_mask, 0, tp5_calc)
            
            tp = tp1 + tp2 + tp3 + tp4 + tp5
            
            # Handle sin_a=0 case if alpha=0 (which is handled in the else block)
            Tav = (ts + tp) / (2 * sin_a**2)
            
            # Recalculate for specific cases where intermediate results might be NaN (like for n=1)
            if np.any(nm_mask):
                # When n=1 (nm=0), the original formula has singularities.
                # Here, we keep the calculated values, and the caller handles the alfa=0 case.
                pass
            
    else: # alfa = 0
        # Tav = 4 * nr / ((nr + 1) * (nr + 1));
        Tav = 4 * nr / ((nr + 1)**2)
        
    return Tav

def soilwat(rdry, nw, kw, SMp, SMC, deleff):
    """
    Translation of MATLAB's soilwat function.
    Calculates wet soil reflectance based on a Poisson process for water films.
    """
    k = np.arange(0, 7) # 0:6
    nk = len(k)
    
    # Mu-parameter of Poisson distribution
    mu = (SMp - 5) / SMC
    
    if np.all(mu <= 0):
        rwet = rdry
    else:
        # Lekner & Dorf (1988) modified soil background reflectance
        # rbac = 1 - (1-rdry) .* (rdry .* tav(90,2.0./nw) / tav(90,2.0) + 1-rdry);
        tav_90_2_nw = _tav(90, 2.0 / nw)
        tav_90_2 = _tav(90, 2.0)
        
        # Need to handle broadcasting for scalar rdry/nw/kw and vector SMp (through mu)
        # However, the MATLAB code structure suggests rdry, nw, kw are [NW, 1] and SMp is [1, NS].
        # In Python, we ensure proper broadcasting: rdry, nw, kw are [N_WL], mu is [N_SMC].
        # rbac will be [N_WL].
        
        rbac = 1 - (1 - rdry) * (rdry * tav_90_2_nw / tav_90_2 + 1 - rdry)
        
        # total reflectance at bottom of water film surface
        # p = 1 - tav(90,nw) ./ nw.^2;
        p = 1 - _tav(90, nw) / nw**2
        
        # reflectance of water film top surface
        # Rw = 1 - tav(40,nw);
        Rw = 1 - _tav(40, nw)
        
        # fractional areas (Poisson distribution)
        # fmul is [nk, NS] if mu is [NS], or [nk] if mu is scalar.
        # np.exp(-mu) * (mu**k / factorial(k))
        
        # Ensure mu is broadcastable with k
        mu_expanded = mu[np.newaxis, :] if mu.ndim == 1 else mu
        k_expanded = k[:, np.newaxis]
        
        # Calculate factorials: 0! to 6!
        factorial_k = np.array([np.math.factorial(i) for i in k])[:, np.newaxis]
        
        # P(k) = mu^k * exp(-mu) / k!
        fmul = np.exp(-mu_expanded) * (mu_expanded**k_expanded / factorial_k)
        
        # two-way transmittance: exp(-2*kw*k*deleff)
        # kw is [N_WL, 1], k is [nk, 1], deleff is scalar. The result is [N_WL, nk]
        kw_expanded = kw[:, np.newaxis]
        
        tw = np.exp(-2 * kw_expanded * deleff * k) # [N_WL, nk]
        
        # Rwet_k: Reflectance of k-thick water film. [N_WL, nk]
        # Rwet_k = Rw + (1-Rw) .* (1-p) .*tw.* rbac ./(1 - p .*tw.* rbac);
        # Ensure proper broadcasting for all terms to [N_WL, nk]
        Rw_exp = Rw[:, np.newaxis]
        p_exp = p[:, np.newaxis]
        rbac_exp = rbac[:, np.newaxis]
        
        numerator = Rw_exp + (1 - Rw_exp) * (1 - p_exp) * tw * rbac_exp
        denominator = 1 - p_exp * tw * rbac_exp
        
        Rwet_k = numerator / denominator
        
        # rwet = rdry * fmul(1) + Rwet_k(:,2:nk)*fmul(2:nk);
        # rwet is [N_WL, NS]
        rdry_exp = rdry[:, np.newaxis]
        
        # Term 1: Dry soil area P(0)
        term1 = rdry_exp * fmul[0, :]
        
        # Term 2: Sum of P(1) to P(nk-1)
        # (N_WL, nk-1) @ (nk-1, NS) -> (N_WL, NS)
        Rwet_k_wet = Rwet_k[:, 1:nk] # P(1) to P(6) reflectance
        fmul_wet = fmul[1:nk, :] # P(1) to P(6) fractional area
        
        # np.matmul requires explicit shape handling. Use broadcasting and summation.
        term2 = np.sum(Rwet_k_wet[:, :, np.newaxis] * fmul_wet[np.newaxis, :, :], axis=1)
        
        rwet = term1 + term2
        
    return rwet


def BSM(soilpar, spec, emp):
    """
    Translation of MATLAB's BSM function.
    Berry Soil Model (BSM) for wet soil reflectance calculation.
    """
    
    # Spectral parameters
    GSV = spec['GSV']  # Global Soil Vectors spectra (nwl * 3)
    kw = spec['Kw']    # water absorption spectrum [nwl]
    nw = spec['nw']    # water refraction index spectrum [nwl]
    
    # Soil parameters
    B = soilpar['BSMBrightness']  # soil brightness
    lat = soilpar['BSMlat']      # spectral shape latitude (deg)
    lon = soilpar['BSMlon']      # spectral shape longitude (deg)
    # soil moisture volume percentage: [0-1] -> [0-100]
    SMp = soilpar['SMC'] * 1E2   # soil moisture volume percentage (5 - 55)
    
    # Empirical parameters
    SMC = emp['SMC']    # soil moisture capacity parameter
    film = emp['film']  # single water film optical thickness
    
    # Convert degrees to radians for sin/cos
    lat_rad = np.deg2rad(lat)
    lon_rad = np.deg2rad(lon)
    
    # Dry soil factors
    f1 = B * np.sin(lat_rad)
    f2 = B * np.cos(lat_rad) * np.sin(lon_rad)
    f3 = B * np.cos(lat_rad) * np.cos(lon_rad)
    
    # Dry soil reflectance: rdry = f1 * GSV(:,1) + f2 * GSV(:,2) + f3 * GSV(:,3)
    # GSV is [N_WL, 3]
    rdry = f1 * GSV[:, 0] + f2 * GSV[:, 1] + f3 * GSV[:, 2]
    
    # Soil moisture effect
    # kw, nw are expected to be [N_WL, 1] for MATLAB's soilwat, so we reshape them if necessary.
    if kw.ndim == 1:
        kw = kw[:, np.newaxis]
    if nw.ndim == 1:
        nw = nw[:, np.newaxis]
        
    rwet = soilwat(rdry, nw, kw, SMp, SMC, film)
    
    return rwet