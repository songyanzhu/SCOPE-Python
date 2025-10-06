import numpy as np

def heatfluxes(ra, rs, Tc, ea, Ta, e_to_q, Ca, Ci, constants, es_fun, s_fun):    
    """
    Translation of src/fluxes/heatfluxes.m
    """
    
    rhoa, cp = constants['rhoa'], constants['cp']

    # Latent heat of vaporization [J kg-1]. Tc is in Celsius.
    lambda_val = (2.501 - 0.002361 * Tc) * 1e6  
    
    ei = es_fun(Tc) # Saturated vapor pressure at Tc
    s = s_fun(ei, Tc) # Slope of saturated vapor pressure curve

    qi = ei * e_to_q # Absolute humidity at surface
    qa = ea * e_to_q # Absolute humidity in air

    # Latent heat flux [W m-2] (Monin-Obukhov/Dalton's Law analogy)
    lE = rhoa / (ra + rs) * lambda_val * (qi - qa)   
    
    # Sensible heat flux [W m-2] (Ohm's Law analogy)
    H = (rhoa * cp) / ra * (Tc - Ta)           
    
    # Vapour pressure at the leaf surface [hPa]
    ec = ea + (ei - ea) * ra / (ra + rs)         
    
    # CO2 concentration at the leaf surface [umol m-3]
    Cc = Ca - (Ca - Ci) * ra / (ra + rs)        
    
    return lE, H, ec, Cc, lambda_val, s