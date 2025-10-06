import numpy as np

def Sint(y, x):
    """
    Simpson-like integration (Trapezoidal Rule applied to bins)
    
    x and y must be vectors (rows, columns) of the same length.
    x must be a monotonically increasing series.
    
    NOTE: The MATLAB code implements the trapezoidal rule, not Simpson's rule.
          The Simpson rule name is a comment error in the original file.
          The function calculates the integral as the sum of areas of trapezoids:
          Integral ≈ Σ [ (y_i + y_{i+1}) / 2 * (x_{i+1} - x_i) ]

    :param y: Function values.
    :param x: Points at which function is evaluated.
    :return: The numerical integral.
    """
    nx = len(x)

    # MATLAB's reshape to ensure column/row vector consistency is handled by NumPy broadcasting,
    # but we'll ensure they are 1D arrays for clear slicing.
    x = np.asarray(x).flatten()
    y = np.asarray(y)
    # If y is 1D, flatten it. If it's 2D (like a row vector in MATLAB), keep shape for matrix multiplication
    if y.ndim == 1:
        y = y.flatten()

    # step = x(2:nx) - x(1:nx-1);  (Differences between adjacent x values)
    step = x[1:nx] - x[0:nx-1]

    # mean = .5 * (y(:,1:nx-1) + y(:,2:nx)); (Mean y value for each segment)
    if y.ndim == 1:
        mean_y = 0.5 * (y[0:nx-1] + y[1:nx])
        # int = mean * step (Dot product for 1D arrays)
        int_val = np.dot(mean_y, step)
    else:
        # Handle 2D y (e.g., [1, nx]) with matrix multiplication (if y is a row-major 2D array)
        # Note: If y is a column vector (nx, 1), the MATLAB operation is different.
        # Assuming y is (n_something, nx) or (1, nx) from the indexing `y(:,1:nx-1)`
        mean_y = 0.5 * (y[:, 0:nx-1] + y[:, 1:nx])
        
        # The MATLAB operation `mean * step` is matrix multiplication: 
        # (n_something, nx-1) * (nx-1, 1) to get (n_something, 1)
        # We assume `step` needs to be a column vector for matrix multiplication in NumPy context
        int_val = np.dot(mean_y, step)
        
    return int_val