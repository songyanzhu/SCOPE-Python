import numpy as np
import os
import warnings

def savebrdfoutput(options, directional, angles, spectral, Output_dir):
    """
    Translation of src/IO/savebrdfoutput.m
    Saves directional output (BRDF, angles, radiances) to ASCII files.
    """

    # --- Setup Directory ---
    directional_dir = os.path.join(Output_dir, 'directional')
    os.makedirs(directional_dir, exist_ok=True)
    
    # --- Prepare Data ---
    
    tts_deg = angles['tts']
    
    # MATLAB: Output_angle = [directional.tto'; directional.psi'];
    Output_angle = np.vstack([directional['tto'], directional['psi']])
    
    # MATLAB: [spectral.wlS directional.refl_]
    Output_refl = np.hstack([spectral['wlS'], directional['refl_']])
    Output_rso = np.hstack([spectral['wlS'], directional['rso_']])
    
    # Files to save and their corresponding data/header
    save_list = [
        (f"refl (SunAngle {tts_deg:.2f} degrees).dat", Output_refl),
        (f"rso (SunAngle {tts_deg:.2f} degrees).dat", Output_rso),
        (f"Angles (SunAngle {tts_deg:.2f} degrees).dat", Output_angle),
    ]

    if options.get('calc_planck', False):
        # MATLAB: [spectral.wlT' directional.Lot_]
        # wlT is 1D array, directional.Lot_ is [n_wlT, n_angles]
        Output_rad = np.hstack([spectral['wlT'].reshape(-1, 1), directional['Lot_']])
        save_list.append((f"Thermal radiances (SunAngle {tts_deg:.2f} degrees).dat", Output_rad))
        
    if options.get('calc_fluor', False):
        # MATLAB: [spectral.wlF' directional.LoF_]
        Output_fluor = np.hstack([spectral['wlF'].reshape(-1, 1), directional['LoF_']])
        save_list.append((f"Fluorescence (SunAngle {tts_deg:.2f} degrees).dat", Output_fluor))

    # --- Save Files (using NumPy for ASCII/TAB separated output) ---
    for filename, data in save_list:
        try:
            np.savetxt(os.path.join(directional_dir, filename), data, 
                       fmt='%.8f', delimiter='\t', header=f"Data for {filename.replace('.dat', '')}", comments='% ')
        except Exception as e:
            warnings.warn(f"Failed to save directional file {filename}: {e}")

    # --- Write read me file ---
    readme_path = os.path.join(directional_dir, 'read me.txt')
    try:
        with open(readme_path, 'w') as fiddirtir:
            fiddirtir.write('The Directional data is written in several files: \r\n')
            fiddirtir.write('\r\n- Angles: contains the observation angles. \r\n')
            fiddirtir.write('   * The 1st row gives the observation zenith angles\r\n')
            fiddirtir.write('   * The 2nd row gives the observation azimuth angles\r\n')
            fiddirtir.write('\r\n- Thermal radiances: contains the radiance at spectral.wlT, emitted plus reflected incoming \r\n')
            fiddirtir.write('   * The 1st column gives the wl values \r\n')
            fiddirtir.write('   * The 2nd column gives the radiances corresponding to the directions given by first column in the Angles file\r\n')
            fiddirtir.write('\r\n- refl and rso: contains the bidirectional distribution functions values, reflectance (based on rso and rdo) and rso. \r\n')
            fiddirtir.write('   * The 1st column gives the wl values corresponding to the BRDF values\r\n')
            fiddirtir.write('   * The 2nd column gives the reflectance values corresponding to the directions given by first column in the Angles file\r\n')
    except IOError as e:
        warnings.warn(f"Failed to write read me.txt: {e}")