import numpy as np

def Monin_Obukhov(constants, meteo, H):
    """
    Calculates the Monin-Obukhov length (L).
    
    :param constants: Structure/dict with rhoa, cp, kappa, g.
    :param meteo: Structure/dict with ustar, Ta.
    :param H: Sensible heat flux (H).
    :return: Monin-Obukhov length L [m].
    """
    # L = -constants.rhoa * constants.cp * meteo.ustar.^3 * ...
    #     (meteo.Ta + 273.15) / (constants.kappa * constants.g * H);
    
    # Use numpy for vectorized operations
    # Note: ustar is expected to be a numpy array for power operation
    # Ta+273.15 is Kelvin
    L = -constants['rhoa'] * constants['cp'] * (meteo['ustar']**3) * \
        (meteo['Ta'] + 273.15) / (constants['kappa'] * constants['g'] * H)
    
    # L(isnan(L)) = -1E6;
    # Replace NaN values with -1E6
    L[np.isnan(L)] = -1E6
    
    return L