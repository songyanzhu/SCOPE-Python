import numpy as np

def satvap(T):
    """
    Calculates the saturated vapor pressure (es) at temperature T [degrees C].
    
    :param T: Temperature [degrees C].
    :return: Saturated vapor pressure es [mbar or hPa].
    """
    
    # constants
    a = 7.5
    b = 237.3  # degrees C

    # calculations
    # es(T) = es(0) * 10^(aT / (b + T))
    # where es(0) = 6.107 mb (hPa)
    es = 6.107 * (10.0**(a * T / (b + T)))
    
    # The MATLAB code had the derivative calculation commented out/set to 0, 
    # but the derivative calculation is present in `slope_satvap.py`.
    
    return es