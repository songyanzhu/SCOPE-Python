import numpy as np

def zo_and_d(soil, canopy, constants):
    """
    Calculates roughness length for momentum (zom) and zero plane displacement (d)
    from vegetation height and LAI using the Verhoef, McNaughton & Jacobs (1997) method.
    
    :param soil: Structure/dict with CSSOIL.
    :param canopy: Structure/dict with CR, CD1, Psicor, LAI, hc.
    :param constants: Structure/dict with kappa.
    :return: A tuple (zom, d) in meters.
    """
    
    # constants
    kappa = constants['kappa']

    # parameters
    CR = canopy['CR']
    CSSOIL = soil['CSSOIL']
    CD1 = canopy['CD1']
    Psicor = canopy['Psicor']
    LAI = canopy['LAI']
    h = canopy['hc']

    # calculations
    # sq = sqrt(CD1 * LAI / 2)
    sq = np.sqrt(CD1 * LAI / 2)
    
    # G1 = max(3.3, (CSSOIL + CR * LAI / 2)^(-0.5))
    G1_base = CSSOIL + CR * LAI / 2
    # Handle zero/negative argument to power if it occurs, though G1_base is usually positive
    G1_power = np.power(G1_base, -0.5, out=np.zeros_like(G1_base), where=G1_base > 0)
    G1 = np.maximum(3.3, G1_power)

    # Zero plane displacement (d) - Eq 12 in Verhoef et al (1997)
    # d = (LAI > 1E-7 & h > 1E-7) * h * (1 - (1 - exp(-sq)) / sq)
    mask = (LAI > 1e-7) & (h > 1e-7)
    
    # Calculate the core term (1 - (1 - exp(-sq)) / sq)
    # Handle the division by sq for non-masked values.
    sq_term = np.zeros_like(sq)
    valid_sq = sq[mask]
    sq_term[mask] = 1 - np.divide(1 - np.exp(-valid_sq), valid_sq, 
                                  out=np.zeros_like(valid_sq), where=valid_sq != 0)
    
    d = h * sq_term * mask # Final d is masked to zero if LAI or h is near zero

    # Roughness length for momentum (zom)
    # zom = (h - d) * exp(-kappa * G1 + Psicor)
    zom = (h - d) * np.exp(-kappa * G1 + Psicor)
    
    return zom, d