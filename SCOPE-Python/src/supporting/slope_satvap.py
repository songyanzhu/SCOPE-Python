import numpy as np

def slope_satvap(T):
    """
    Calculates the saturated vapor pressure (es) and its derivative (s) 
    with respect to temperature T [degrees C].
    
    :param T: Temperature [degrees C].
    :return: A tuple (es, s) where es is saturated vapor pressure [mbar or hPa] 
             and s is the derivative [mbar/C or hPa/C].
    """
    
    # constants
    a = 7.5
    b = 237.3  # degrees C
    log10 = 2.3026 # ln(10)

    # calculations
    # es(T) = 6.107 * 10^(aT / (b + T))
    es = 6.107 * (10.0**(a * T / (b + T)))
    
    # s(T) = es * ln(10) * a * b / (b + T)^2
    s = es * log10 * a * b / (b + T)**2
    
    return es, s