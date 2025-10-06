import numpy as np

def count_k(nvars, v, vmax, id_start):
    """
    Increments a multi-digit counter vector, resetting digits at their max.
    
    :param nvars: Number of digits (length of v and vmax).
    :param v: Current vector of digits (1-based count in MATLAB).
    :param vmax: Maximum values for each digit.
    :param id_start: Starting digit index (1-based in MATLAB).
    :return: New vector of digits (vnew).
    """
    
    # Convert inputs to 0-based indexing for Python
    # Make copies to avoid modifying inputs
    v_py = np.asarray(v, dtype=int).copy()
    vmax_py = np.asarray(vmax, dtype=int).copy()
    i = id_start - 1 # 0-based index
    
    # The logic is based on 1-based indices in MATLAB, 
    # so we must convert to 0-based for list/array access.
    # v(i) == vmax(i)
    while v_py[i] == vmax_py[i]:
        # v(i) = 1;
        v_py[i] = 1
        
        # i = rem(i, nvars) + 1; (move to the next index, wrapping around)
        i = (i + 1) % nvars
        
    # v(i) = rem(v(i), vmax(i)) + 1; 
    # This is equivalent to v(i) = v(i) + 1, since the loop exits 
    # when v(i) < vmax(i). If v(i) = vmax(i), it would have been reset to 1.
    # The MATLAB modulo logic for incrementing (when not max) is:
    # If v(i) = 1 and vmax(i) = 5, rem(1, 5) + 1 = 2 (correct)
    # If v(i) = 4 and vmax(i) = 5, rem(4, 5) + 1 = 5 (correct)
    # If v(i) = 5 and vmax(i) = 5, rem(5, 5) + 1 = 1 (this case should not be reached)
    
    v_py[i] = v_py[i] + 1
    
    return v_py