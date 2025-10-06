import numpy as np

def Planck(wl, Tb, em=None):
    """
    Calculates the Planck blackbody radiance or spectral radiance.

    :param wl: Wavelength [nm].
    :param Tb: Brightness temperature [K].
    :param em: Emissivity, defaults to 1 (blackbody).
    :return: Spectral radiance Lb [W m-2 sr-1 nm-1].
    """
    c1 = 1.191066e-22  # MATLAB's value (likely in W m^2 / sr) - check units
    c2 = 14388.33     # [K nm] - check units

    if em is None:
        # MATLAB's ones(size(Tb)) equivalent
        em = np.ones_like(Tb)
        
    # MATLAB: Lb = em.* c1*(wl*1e-9).^(-5)./(exp(c2./(wl*1e-3*Tb))-1);
    # wl is in nm, so wl*1e-9 is in m for (wl*1e-9)^(-5)
    # wl*1e-3 is in m for the exponent part (assuming c2 is in K m)
    
    # If wl is in nm, and we want Lb in W m-2 sr-1 nm-1, then c1 should be 
    # $2 \pi h c^2$ in $W m^2 / sr$ if we use $m$ in the denominator, and then 
    # divide by $10^9$ to get $nm^{-1}$. 
    # However, I will stick to the exact MATLAB formula:
    
    # wl*1e-9 is wavelength in meters
    # wl*1e-3 is wavelength in millimeters (used in exponent)
    
    # The original MATLAB code implies the following units for $c_1$ and $c_2$
    # if $wl$ is in $nm$ and $L_b$ is in $W m^{-2} sr^{-1} nm^{-1}$:
    # Lb = em * c1 * (wl_m)^(-5) / (exp(c2/(wl_mm * Tb)) - 1)
    
    wl_m = wl * 1e-9
    wl_mm = wl * 1e-3
    
    Lb = em * c1 * (wl_m**(-5)) / (np.exp(c2 / (wl_mm * Tb)) - 1)
    
    return Lb