# Import the function it calls
from .ephoton import ephoton
import numpy as np

def e2phot(lambda_m, E, constants):
    """
    Calculates the number of moles of photons corresponding to E Joules of energy
    of wavelength lambda [m].
    
    :param lambda_m: Wavelength [m].
    :param E: Energy [J].
    :param constants: Structure/dict with h, c, A.
    :return: molphotons [mol].
    """
    
    # e = ephoton(lambda, constants); (Energy of a single photon)
    e = ephoton(lambda_m, constants)
    
    # photons = E / e;
    photons = E / e
    
    # molphotons = photons / constants.A; (constants.A is Avogadro's number)
    molphotons = photons / constants['A']
    
    return molphotons