import numpy as np

def meanleaf(canopy, F, choice, Ps=None):
    """
    Calculates the layer average and the canopy average of leaf properties.

    :param canopy: Structure/dict with nlayers, nlincl, nlazi, lidf.
    :param F: Input matrix [nlincl, nlazi, nlayers].
    :param choice: Integration method ('angles', 'layers', or 'angles_and_layers').
    :param Ps: Fraction sunlit per layer [nlayers] (optional, required for 'layers'/'angles_and_layers').
    :return: The averaged quantity Fout.
    """
    
    # Extract dimensions and LIDF
    nl = canopy['nlayers']
    nli = canopy['nlincl']
    nlazi = canopy['nlazi']
    lidf = np.asarray(canopy['lidf']) # Ensure lidf is a numpy array
    
    F = np.asarray(F) # Ensure F is a numpy array
    
    # Initialize Fout (it's often reshaped/reduced later)
    # The initialization of Fout in MATLAB is `Fout = zeros(nli, nlazi,nl);` 
    # but Fout is reassigned fully within the switch cases.
    Fout = None

    if choice == 'angles':
        # Integration over leaf angles
        # lidf is assumed to be a vector of size nli
        
        # Fout(j,:,:) = F(j,:,:)*lidf(j);
        # Use broadcasting: lidf needs to be (nli, 1, 1)
        lidf_reshaped = lidf.reshape(nli, 1, 1)
        F_weighted = F * lidf_reshaped
        
        # Fout = sum(sum(Fout))/nlazi; [1,1,nl] - sum over nli and nlazi
        Fout = np.sum(F_weighted, axis=(0, 1)) / nlazi
        
        # Fout = permute(Fout,[3 1 2]); [nl] 
        # Fout is already (nl,) after the sum reduction (np.sum axes 0, 1)
        # and doesn't need permutation if it's already 1D.
        
    elif choice == 'layers':
        # Integration over layers only
        # Fout = Ps'*F/nl; (Dot product of row vector Ps' and matrix F)
        if Ps is None:
            raise ValueError("Ps is required for 'layers' choice.")
        Ps = np.asarray(Ps).flatten()
        Fout = np.dot(Ps, F) / nl
        
    elif choice == 'angles_and_layers':
        # Integration over both leaf angles and layers
        
        # Weight by LIDF (leaf angles)
        lidf_reshaped = lidf.reshape(nli, 1, 1)
        F_weighted_angles = F * lidf_reshaped
        
        # Weight by sunlit fraction Ps (layers)
        if Ps is None:
            raise ValueError("Ps is required for 'angles_and_layers' choice.")
        Ps = np.asarray(Ps).flatten() # ensure Ps is 1D
        
        # Fout(:,:,j) = Fout(:,:,j)*Ps(j);
        # Use broadcasting: Ps needs to be (1, 1, nl)
        Ps_reshaped = Ps.reshape(1, 1, nl)
        F_weighted_both = F_weighted_angles * Ps_reshaped
        
        # Fout = sum(sum(sum(Fout))) / nlazi / nl; [1]
        Fout = np.sum(F_weighted_both) / nlazi / nl
        
    else:
        raise ValueError(f"Invalid choice for meanleaf: {choice}")
        
    return Fout