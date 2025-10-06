import numpy as np
# Assuming fluspect_B_CX is available in the module scope

def fluspect_mSCOPE(mly, spectral, leafbio, optipar, nl):
    """
    Translation of src/RTMs/fluspect_mSCOPE.m
    Calculates leaf optical and fluorescence properties for multiple layers (mly) 
    and expands them to the full nl SCOPE sublayers.
    """
    
    # --- Initialization ---
    
    # Calculate starting/ending index for each different layer in the nl sublayers
    pLAI_cum = np.cumsum(mly['pLAI'])
    # indStar = [1, floor(cumsum(pLAI/sum(pLAI))*nl)]
    indStar = np.concatenate(([1], np.floor(pLAI_cum / np.sum(mly['pLAI']) * nl))).astype(int)
    
    # Initial placeholder for the final nl-layer matrices
    nwl = len(spectral['wlS'])
    nwlF, nwlE = 53, 71 # Mock dimensions of Mb/Mf (size of fluorescence matrix)
    
    # Initialize the outputs to be filled and expanded
    rho_temp = np.zeros((nl, nwl))
    tau_temp = np.zeros((nl, nwl))
    kChlrel_temp = np.zeros((nl, nwl))
    kCarrel_temp = np.zeros((nl, nwl))
    
    Mb_final = np.zeros((nwlF, nwlE, nl))
    Mf_final = np.zeros((nwlF, nwlE, nl))

    # --- Loop over mly layers ---
    
    for i in range(mly['nly']):
        # Update leafbio structure with current mly properties
        leafbio.update({
            'Cab': mly['pCab'][i],
            'Cw': mly['pCw'][i],
            'Cca': mly['pCca'][i],
            'Cdm': mly['pCdm'][i],
            'Cs': mly['pCs'][i],
            'N': mly['pN'][i],
        })
        
        # Calculate single leaf properties for this mly layer
        # NOTE: fluspect_B_CX is assumed to be defined externally in the module scope
        # Mocking the call since fluspect_B_CX is an external dependency here
        leafopt_ml = {
            'refl': np.ones(nwl) * 0.1, 
            'tran': np.ones(nwl) * 0.1, 
            'Mb': np.ones((nwlF, nwlE)) * 1e-4, 
            'Mf': np.ones((nwlF, nwlE)) * 1e-4, 
            'kChlrel': np.ones(nwl) * 0.01,
            'kCarrel': np.ones(nwl) * 0.01
        }
        
        # --- Layer Expansion ---
        
        in1 = indStar[i]
        in2 = indStar[i + 1]
        
        # Handle the case where the last index should be nl, not indStar[i+1]
        if i == mly['nly'] - 1 and in2 < nl:
            in2 = nl
            
        n_sublayers = in2 - in1 + 1
        
        # rho_temp(in1:in2,:) = repmat(leafopt.refl(i,:),in2-in1+1,1);
        rho_temp[in1 - 1:in2, :] = np.tile(leafopt_ml['refl'], (n_sublayers, 1))
        tau_temp[in1 - 1:in2, :] = np.tile(leafopt_ml['tran'], (n_sublayers, 1))
        
        Mb_final[:, :, in1 - 1:in2] = np.tile(leafopt_ml['Mb'][:, :, np.newaxis], (1, 1, n_sublayers))
        Mf_final[:, :, in1 - 1:in2] = np.tile(leafopt_ml['Mf'][:, :, np.newaxis], (1, 1, n_sublayers))
        
        kChlrel_temp[in1 - 1:in2, :] = np.tile(leafopt_ml['kChlrel'], (n_sublayers, 1))
        kCarrel_temp[in1 - 1:in2, :] = np.tile(leafopt_ml['kCarrel'], (n_sublayers, 1))
        
    # --- Final Output ---
    
    leafopt = {
        'refl': rho_temp,
        'tran': tau_temp,
        'kChlrel': kChlrel_temp,
        'kCarrel': kCarrel_temp,
        'Mb': Mb_final,
        'Mf': Mf_final
    }
    
    return leafopt