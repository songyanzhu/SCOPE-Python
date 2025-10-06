import numpy as np

def dcum(a, b, theta):
    """
    Subroutine for Cumulative Leaf Inclination Density Function (dcum).
    
    :param a: Parameter a.
    :param b: Parameter b.
    :param theta: Leaf angle [degrees].
    :return: Cumulative Leaf Inclination Density Function F.
    """
    rd = np.pi / 180  # Geometrical constant (deg to rad)
    
    if a > 1:
        # F = 1 - cos(theta * rd)
        F = 1 - np.cos(theta * rd)
    else:
        # Note: This is an iterative solver (Newton-Raphson-like)
        eps = 1e-8
        delx = 1.0  # Initial large value
        
        x = 2 * rd * theta
        theta2 = x.copy() # theta2 = x
        
        # Vectorized implementation of the while loop (if inputs are arrays)
        # Using a fixed max iteration for robustness if inputs are arrays and convergence is mixed
        # For a pure scalar input, a standard while loop would be fine.
        
        # Assuming scalar inputs for simplicity matching the MATLAB context
        max_iter = 100
        for _ in range(max_iter):
            y = a * np.sin(x) + 0.5 * b * np.sin(2 * x)
            dx = 0.5 * (y - x + theta2)
            
            x = x + dx
            delx = np.abs(dx)
            
            if delx < eps:
                break
            
        F = (2 * y + theta2) / np.pi 
        
    return F

def leafangles(a, b):
    """
    Calculates the Leaf Inclination Distribution Function (LIDF) vector.

    :param a: Parameter a for the distribution.
    :param b: Parameter b for the distribution.
    :return: LIDF vector (13 elements).
    """
    
    F = np.zeros(13)
    
    # i=1:8 (0-based: 0 to 7)
    for i in range(8):
        theta = (i + 1) * 10  # 10, 20, ..., 80 degrees
        F[i] = dcum(a, b, theta)
        
    # i=9:12 (0-based: 8 to 11)
    for i in range(8, 12):
        theta = 80 + (i - 7) * 2  # 82, 84, 86, 88 degrees
        F[i] = dcum(a, b, theta)
        
    # i=13:13 (0-based: 12)
    F[12] = 1  # For 90 degrees
    
    lidf = np.zeros(13)
    
    # i=13:-1:2 (0-based: 12 down to 1)
    for i in range(12, 0, -1):
        # lidf(i) = F(i) - F(i-1); (MATLAB 1-based to Python 0-based)
        lidf[i] = F[i] - F[i-1]
        
    # lidf(1) = F(1); (0-based index 0)
    lidf[0] = F[0]
    
    return lidf