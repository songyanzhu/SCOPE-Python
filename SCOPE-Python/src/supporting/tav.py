import numpy as np

def tav(alfa, nr):
    """
    Calculates the average transmittance (Tav) of a dielectric interface
    for unpolarized light, likely using an approximation derived from Fresnel's equations.
    (This function implements the cubic equation method by Lekner & Dorf).

    :param alfa: Incident angle [degrees].
    :param nr: Refractive index.
    :return: Average transmittance Tav.
    """
    
    n2 = nr**2
    np_val = n2 + 1
    nm = n2 - 1

    alfa_rad = np.deg2rad(alfa)
    sin_a = np.sin(alfa_rad)
    
    if alfa != 0:
        # B2 = sin_a^2 - np/2
        B2 = sin_a**2 - np_val / 2
        
        # k = -((n2 - 1)^2) / 4
        k = -(nm**2) / 4
        
        # B1 = (alfa != 90) * sqrt(B2^2 + k)
        # Using a mask for alfa=90 is simpler in NumPy
        B1_base = np.sqrt(B2**2 + k)
        if alfa == 90:
            B1 = 0.0 # B1 = 0 if alfa = 90
        else:
            B1 = B1_base

        # b = B1 - B2
        b = B1 - B2
        
        # a = +((nr + 1)^2) / 2
        a = ((nr + 1)**2) / 2
        
        b3 = b**3
        a3 = a**3
        
        # ts (Transmittance S-polarized) approximation
        # ts = (k^2 / (6*b3) + k / b - b / 2) - ...
        #      (k^2 / (6*a3) + k / a - a / 2)
        ts_term1 = np.divide(k**2, 6 * b3, out=np.zeros_like(b3), where=b3 != 0) + np.divide(k, b, out=np.zeros_like(b), where=b != 0) - b / 2
        ts_term2 = np.divide(k**2, 6 * a3, out=np.zeros_like(a3), where=a3 != 0) + np.divide(k, a, out=np.zeros_like(a), where=a != 0) - a / 2
        ts = ts_term1 - ts_term2
        
        # tp (Transmittance P-polarized) approximation terms (tp1 to tp5)
        
        # tp1 = -2 * n2 * (b - a) / (np^2)
        tp1 = -2 * n2 * (b - a) / (np_val**2)
        
        # tp2 = -2 * n2 * np * log(b / a) / (nm^2)
        # Handle log(0) and div by zero
        log_ba = np.log(np.divide(b, a, out=np.zeros_like(b), where=a != 0))
        tp2 = -2 * n2 * np_val * log_ba / (nm**2)
        
        # tp3 = n2 * (1 / b - 1 / a) / 2
        tp3 = n2 * (np.divide(1, b, out=np.zeros_like(b), where=b != 0) - np.divide(1, a, out=np.zeros_like(a), where=a != 0)) / 2
        
        # The denominators in tp4 and tp5 are nm^2 and np^3.
        term_denom = 2 * np_val * b - nm**2
        
        # tp4 = 16 * n2^2 * (n2^2 + 1) * log((2*np*b - nm^2) / (2*np*a - nm^2)) / (np^3 * nm^2)
        tp4_num = np.log(np.divide(2 * np_val * b - nm**2, 2 * np_val * a - nm**2, 
                                 out=np.zeros_like(b), where=(2 * np_val * a - nm**2) != 0))
        tp4 = 16 * n2**2 * (n2**2 + 1) * tp4_num / (np_val**3 * nm**2)
        
        # tp5 = 16 * n2^2 * n2 * (1/(2*np*b - nm^2) - 1/(2*np*a - nm^2)) / (np^3)
        tp5_term = np.divide(1, 2 * np_val * b - nm**2, out=np.zeros_like(b), where=(2 * np_val * b - nm**2) != 0) - \
                   np.divide(1, 2 * np_val * a - nm**2, out=np.zeros_like(a), where=(2 * np_val * a - nm**2) != 0)
        tp5 = 16 * n2**2 * n2 * tp5_term / (np_val**3)
        
        tp = tp1 + tp2 + tp3 + tp4 + tp5
        
        # Tav = (ts + tp) / (2 * sin_a^2)
        Tav = np.divide(ts + tp, 2 * sin_a**2, 
                        out=np.zeros_like(sin_a), where=sin_a != 0)

    else:
        # Normal incidence (alfa = 0)
        # Tav = 4 * nr / ((nr + 1) * (nr + 1))
        Tav = 4 * nr / ((nr + 1)**2)
        
    return Tav