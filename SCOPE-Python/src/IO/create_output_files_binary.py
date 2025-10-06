import numpy as np
import os
import shutil
import time
import warnings

def create_output_files_binary(parameter_file, F, path_of_code, input_path, spectral, options):
    """
    Translation of src/IO/create_output_files_binary.m
    Creates the output directory and opens binary files for writing.
    """
    
    # MATLAB: string = clock; simulation_name = char(F(1).FileName);
    sim_time = time.localtime()
    simulation_name = F[0]['FileName'] 
    
    # MATLAB: outdir_name = sprintf('%s_%4.0f-%02.0f-%02.0f-%02.0f%02.0f', simulation_name, string(1:5));
    outdir_name = f"{simulation_name}_{sim_time.tm_year:04d}-{sim_time.tm_mon:02d}-{sim_time.tm_mday:02d}-{sim_time.tm_hour:02d}{sim_time.tm_min:02d}"
    
    # MATLAB: Output_dir = [fullfile('output', outdir_name) filesep];
    Output_dir = os.path.join('output', outdir_name) + os.sep

    # %% Create Output dir
    # MATLAB: warning('off','MATLAB:DELETE:FileNotFound') and mkdir logic
    try:
        os.makedirs(Output_dir, exist_ok=True)
        os.makedirs(os.path.join(Output_dir, 'Parameters'), exist_ok=True)
    except Exception as e:
        raise IOError(f"Failed to create output directory {Output_dir}: {e}")

    # %% Log File (Parameter files copy)
    for p_file in parameter_file[0]:
        # MATLAB: copy_name = [strrep(parameter_file{1}{i}, '.csv', '') '_' outdir_name '.csv'];
        base_name = p_file.replace('.csv', '')
        copy_name = f"{base_name}_{outdir_name}.csv"
        
        src = os.path.join(input_path, p_file)
        dst = os.path.join(Output_dir, 'Parameters', copy_name)
        
        try:
            shutil.copy2(src, dst)
        except FileNotFoundError:
            warnings.warn(f"Input parameter file not found: {src}")

    # Write SCOPE version/path
    fidpath = os.path.join(Output_dir, 'Parameters', 'SCOPEversion.txt')
    try:
        with open(fidpath, 'w') as f_log:
            f_log.write(path_of_code)
    except IOError as e:
        warnings.warn(f"Failed to write SCOPEversion.txt: {e}")

    # %% Filenames
    fnames = {}
    fnames['pars_file'] = os.path.join(Output_dir, 'pars_and_input_short.bin')
    fnames['apar_file'] = os.path.join(Output_dir, 'aPAR.bin')
    fnames['veg_file'] = os.path.join(Output_dir, 'vegetation.bin')
    fnames['flu_file'] = os.path.join(Output_dir, 'fluxes.bin')
    fnames['rad_file'] = os.path.join(Output_dir, 'radiation.bin')
    
    if options['calc_fluor']:
        fnames['fluor_file'] = os.path.join(Output_dir, 'fluorescence_scalars.bin')
        fnames['fluor_spectrum_file'] = os.path.join(Output_dir, 'fluorescence.bin')
        fnames['sigmaF_file'] = os.path.join(Output_dir, 'sigmaF.bin')
        fnames['fhemis_file'] = os.path.join(Output_dir, 'fluorescence_hemis.bin')
        fnames['fRC_file'] = os.path.join(Output_dir, 'fluorescence_ReabsCorr.bin')
        fnames['fRCL_file'] = os.path.join(Output_dir, 'fluorescence_AllLeaves.bin')
        fnames['Lo2_file'] = os.path.join(Output_dir, 'Lo_spectrum_inclF.bin')
        fnames['rapp_file'] = os.path.join(Output_dir, 'apparent_reflectance.bin')

    if options['save_spectral']:
        fnames['r_file'] = os.path.join(Output_dir, 'reflectance.bin')
        fnames['rsd_file'] = os.path.join(Output_dir, 'rsd.bin')
        fnames['rdd_file'] = os.path.join(Output_dir, 'rdd.bin')
        fnames['rso_file'] = os.path.join(Output_dir, 'rso.bin')
        fnames['rdo_file'] = os.path.join(Output_dir, 'rdo.bin')
        fnames['Eout_file'] = os.path.join(Output_dir, 'Eout_spectrum.bin')
        fnames['Lo_file'] = os.path.join(Output_dir, 'Lo_spectrum.bin')
        fnames['Esun_file'] = os.path.join(Output_dir, 'Esun.bin')
        fnames['Esky_file'] = os.path.join(Output_dir, 'Esky.bin')
        
    fnames['resist_file'] = os.path.join(Output_dir, 'resistances.bin')

    # %% Open files for writing (binary mode 'wb')
    f = {}
    for name, path in fnames.items():
        try:
            f[name] = open(path, 'wb')
        except IOError as e:
            raise IOError(f"Failed to open binary file {path}: {e}")

    # %% write wl
    # MATLAB's save(..., '-ascii') is emulated by numpy.savetxt
    wlS = spectral['wlS']
    wlF = spectral['wlF']

    try:
        np.savetxt(os.path.join(Output_dir, 'wlS.txt'), wlS, fmt='%.1f')
        np.savetxt(os.path.join(Output_dir, 'wlF.txt'), wlF, fmt='%.1f')
    except Exception as e:
        warnings.warn(f"Failed to save wavelength files (wlS/wlF.txt): {e}")

    return Output_dir, f, fnames