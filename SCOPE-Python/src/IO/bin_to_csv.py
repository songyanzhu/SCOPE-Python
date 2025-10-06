import numpy as np
import os
import struct
import warnings

# --- Helper Function (Private to Module) ---

def _write_output(header, units, bin_path, f_n_col, ns, not_header=False):
    """
    Translation of MATLAB's internal write_output function.
    Reads binary data, reshapes it, and writes it to a CSV file.
    """
    
    n_csv = bin_path.replace('.bin', '.csv')

    # Open CSV file for writing ('w')
    try:
        f_csv = open(n_csv, 'w')
    except IOError as e:
        warnings.warn(f"Could not open CSV file {n_csv}: {e}")
        return

    # 1. Write Header
    header_str = ','.join(header) + '\n'
    if not_header:
        header_str = '#' + header_str
    else:
        # Check if number of headers matches column count
        if len(header) != f_n_col:
            raise ValueError(f"Less headers ({len(header)}) than columns ({f_n_col}) in {bin_path}")
            
    f_csv.write(header_str)
    
    # 2. Write Units
    f_csv.write('#' + ','.join(units) + '\n')

    # 3. Read Binary Data
    try:
        with open(bin_path, 'rb') as f_bin:
            # Read all data as 'double' (8 bytes)
            data = f_bin.read()
            # Calculate total number of doubles
            num_doubles = len(data) // 8
            
            # Unpack all doubles
            out = struct.unpack(f'<{num_doubles}d', data)
            out = np.array(out)
            
    except IOError as e:
        warnings.warn(f"Could not read binary file {bin_path}: {e}")
        f_csv.close()
        return

    # 4. Reshape and Write Data Rows
    # MATLAB: out_2d = reshape(out, f_n_col, ns)'
    try:
        out_2d = out.reshape(f_n_col, ns).T
    except ValueError:
        warnings.warn(f"Binary data size {out.size} does not match expected shape ({f_n_col}, {ns}) for {bin_path}. Skipping data write.")
        f_csv.close()
        return

    # MATLAB: dlmwrite/fprintf loop
    for k in range(ns):
        # Format the row, ensuring integer precision as in MATLAB's dlmwrite/fprintf loop ('%d')
        # We use NumPy's tolist() to make it easier to format
        row_list = out_2d[k, :].tolist()
        
        # Format: '%d,' for all but last element, '%d\n' for last element
        formatted_row = ','.join([f"{x:.10g}" for x in row_list[:-1]])
        if formatted_row:
             formatted_row += ','
        formatted_row += f"{row_list[-1]:.10g}\n"
        
        f_csv.write(formatted_row)
        
    f_csv.close()


# --- Main Function ---

def bin_to_csv(fnames, V, vmax, n_col, ns):
    """
    Translation of src/IO/bin_to_csv.m
    Converts binary output files to CSV format.
    """
    
    # MATLAB's structfun(@delete, fnames) is deferred until the end.
    files_to_delete = []

    # Helper function to call _write_output for a set of files and store for deletion
    def process_file(name, header, units, n_col_key, delete_bin=True, not_header=False):
        if name in fnames:
            path = fnames[name]
            _write_output(header, units, path, n_col[n_col_key], ns, not_header)
            if delete_bin:
                files_to_delete.append(path)

    # %% pars
    vmax_indices = np.where(vmax > 1)[0]
    if vmax_indices.size > 0:
        pars_names = ['n_pars'] + [V[i]['Name'] for i in vmax_indices]
        process_file('pars_file', pars_names, [''] * len(pars_names), 'pars', delete_bin=False)

    # %% aPAR
    apar_names = [
        'simulation_number', 'year', 'DoY', 'iPAR', 'iPARE', 'LAIsunlit', 'LAIshaded',
        'aPARtot', 'aPARsun', 'aPARsha', 'aPARCabtot', 'aPARCabsun', 'aPARCabsha',
        'aPARCartot', 'aPARCarsun', 'aPARCarsha', 'aPARtotE', 'aPARsunE', 'aPARshaE',
        'aPARCabtotE', 'aPARCabsunE', 'aPARCabshaE', 'aPARCartotE', 'aPARCarsunE', 'aPARCarshaE',
    ]
    apar_units = [
        '', '', '', 'umol m-2 s-1','W m-2', 'm2m-2','m2m-2',
        'umol m-2 s-1', 'umol m-2 s-1', 'umol m-2 s-1', 'umol m-2 s-1', 'umol m-2 s-1', 'umol m-2 s-1',
        'umol m-2 s-1', 'umol m-2 s-1', 'umol m-2 s-1', 'W m-2','W m-2','W m-2',
        'W m-2','W m-2','W m-2', 'W m-2','W m-2','W m-2'
    ]
    process_file('apar_file', apar_names, apar_units, 'apar')

    # %% veg
    veg_names = ['simulation_number', 'year', 'DoY', 'Photosynthesis', 'Electron_transport', 'NPQ_energy', 'NPQ_photon', 'canopy_level_FQE','LST','emis', 'GPP']
    veg_units = ['', '', '', 'umol CO2 m-2 s-1', 'umol m-2 s-1', 'W m-2', 'umol m-2 s-1', 'umol photons (umol photons)-1', 'K','', 'umol CO2 m-2 s-1']
    process_file('veg_file', veg_names, veg_units, 'veg')

    # %% flu
    flu_names = [
        'simulation_number','nu_iterations', 'year','DoY',
        'Rnctot','lEctot','Hctot','Actot','Tcave', 
        'Rnstot','lEstot','Hstot','Gtot','Tsave',
        'Rntot','lEtot','Htot'
    ]
    flu_units = [
        '', '', '', '', 
        'W m-2','W m-2','W m-2','umol m-2 s-1','C',
        'W m-2','W m-2','W m-2','W m-2','C',
        'W m-2',' W m-2','W m-2'
    ]
    process_file('flu_file', flu_names, flu_units, 'flu')

    # %% rad
    rad_names = ['simulation_number','year','DoY','ShortIn','LongIn','HemisOutShort','HemisOutLong','Lo','Lot','Lote']
    rad_units = ['','','','W m-2','W m-2','W m-2','W m-2','W m-2 sr-1','W m-2 sr-1','W m-2 sr-1']
    process_file('rad_file', rad_names, rad_units, 'rad')

    # %% fluor
    if 'fluor_file' in fnames:
        fluor_names = ['F_1stpeak', 'wl_1stpeak', 'F_2ndpeak', 'wl_2ndpeak', 'F684', 'F761', 'LFtot', 'EFtot', 'EFtot_RC']
        fluor_units = ['W m-2 um-1 sr-1','nm','W m-2 um-1 sr-1','nm','W m-2 um-1 sr-1','W m-2 um-1 sr-1','W m-2 sr-1','W m-2','W m-2']
        process_file('fluor_file', fluor_names, fluor_units, 'fluor')

        process_file('fluor_spectrum_file', ['fluorescence_spectrum 640:1:850 nm'], ['W m-2 um-1 sr-1'], 'fluor_spectrum', not_header=True)
        process_file('sigmaF_file', ['escape probability 640:1:850 nm'], [''], 'sigmaF', not_header=True)
        process_file('fhemis_file', ['fluorescence_spectrum 640:1:850 nm hemispherically integrated'], ['W m-2 um-1'], 'fhemis', not_header=True)
        process_file('fRC_file', ['fluorescence_spectrum 640:1:850 nm reabsorption corrected'], ['W m-2 um-1'], 'fRC', not_header=True)
        process_file('fRCL_file', ['fluorescence_spectrum 640:1:850 nm emission by all leaves'], ['W m-2 um-1'], 'fRCL', not_header=True)
        process_file('Lo2_file', ['upwelling radiance including fluorescence'], ['W m-2 um-1 sr-1'], 'Lo2', not_header=True)
        process_file('rapp_file', ['apparent reflectance'], [''], 'rapp', not_header=True)

    # %% reflectance
    if 'r_file' in fnames:
        process_file('r_file', ['reflectance'], ['pi*upwelling radiance/irradiance'], 'r', not_header=True)
        process_file('rsd_file', ['rsd'], ['directional-hemispherical reflectance factor'], 'rsd', not_header=True)
        process_file('rdd_file', ['rdd'], ['bi-hemispherical reflectance factor'], 'rdd', not_header=True)
        process_file('rso_file', ['rso'], ['bi-directional reflectance factor'], 'rso', not_header=True)
        process_file('rdo_file', ['rdo'], ['hemispherical-directional reflectance factor'], 'rdo', not_header=True)

        # %% radiance
        process_file('Eout_file', ['hemispherically integrated upwelling radiance'], ['W m-2 um-1'], 'Eout', not_header=True)
        process_file('Lo_file', ['upwelling radiance excluding fluorescence'], ['W m-2 um-1 sr-1'], 'Lo', not_header=True)
        process_file('Esun_file', ['direct solar irradiance'], ['W m-2 um-1'], 'Esun', not_header=True)
        process_file('Esky_file', ['diffuse solar irradiance'], ['W m-2 um-1'], 'Esky', not_header=True)

    # %% resistances
    process_file('resist_file', ['aerodynamicresistance(ra)', 'raforsoil', 'rss','ustar'], ['s m-1','s m-1','s m-1','m s-1'], 'resist', not_header=True)

    # MATLAB: fclose('all') is implicit here as all files are closed in _write_output.

    # %% deleting .bin
    for path in files_to_delete:
        try:
            os.remove(path)
        except OSError as e:
            warnings.warn(f"Could not delete binary file {path}: {e}")