import numpy as np

def ephoton(lambda_m, constants):
    """
    Calculates the energy content (J) of 1 photon of wavelength lambda [m].
    
    :param lambda_m: Wavelength [m].
    :param constants: Structure/dict with h (Planck's constant) and c (speed of light).
    :return: Energy of 1 photon E [J].
    """
    
    h = constants['h']  # Planck's constant [J s]
    c = constants['c']  # Speed of light [m s-1]
    
    # E = h * c / lambda;
    E = h * c / lambda_m
    
    return E