import numpy as np
import os
import sys
import math
from datetime import datetime

# --- Import previously translated functions ---
# Assuming 'src/fluxes' is in the Python path or the current directory
try:
    # Adjust path to find modules if running from the root
    sys.path.append(os.path.join(os.path.dirname(__file__), 'src', 'fluxes'))
    from src.fluxes.biochemical import biochemical
    from src.fluxes.biochemical_md12 import biochemical_MD12
    from src.fluxes.heatfluxes import heatfluxes
    from src.fluxes.resistances import resistances
    from src.fluxes.ebal import ebal
except ImportError as e:
    print(f"Warning: Could not import internal SCOPE modules. Using placeholders. Error: {e}")
    # Define local placeholders if imports fail
    def biochemical(*args): return {'A': 0.0, 'Ag': 0.0, 'Ja': 0.0, 'eta': 0.0, 'Kn': 0.0, 'Phi_N': 0.0}
    def biochemical_MD12(*args): return {'A': 0.0, 'Ag': 0.0, 'Ja': 0.0, 'eta': 0.0, 'Kn': 0.0, 'Phi_N': 0.0}
    def ebal(*args): 
        # Mocking complex EBAL output
        nl = args[5]['nlayers']
        return {'counter': 1}, {'canopyemis': 0.95, 'Lot': 1.0, 'Lote': 1.0, 'Eoutte': 1.0, 'LoF_': 0.1, 'Esun_': np.zeros(1800), 'Esky_': np.zeros(1800), 'refl': np.zeros(1800), 'Lo_': np.zeros(1800), 'Lot_': np.zeros(1800)}, {'Tcu': np.zeros(nl), 'Tch': np.zeros(nl), 'Tsu': 295, 'Tsh': 293}, {'Tsold': np.zeros((12, 2))}, {}, {}, {}, {'rss': 50}, args[4]

# --- Placeholder Functions (for external MATLAB dependencies) ---

def define_constants():
    """Mocks define_constants.m"""
    return {
        'MH2O': 18.01528, 'Mair': 28.96, 'R': 8.314, 'rhoa': 1.2047, 'cp': 1004, 
        'sigmaSB': 5.67e-8, 'kappa': 0.41, 'A': 6.02214086e23 # Avogadro
    }

def define_bands():
    """Mocks define_bands.m"""
    # Simplified spectral definition
    wl = np.arange(400, 2501) # 400nm to 2500nm
    wlT = np.arange(10000, 15000) # Thermal
    wlF = np.arange(640, 851) # Fluorescence
    return {
        'wl': wl, 'wlF': wlF, 'wlT': wlT, 
        'IwlP': np.where((wl >= 400) & (wl <= 700))[0],
        'IwlT': np.where((wl >= 700) & (wl <= 2500))[0],
        'SCOPEspec': None # Placeholder for external spectral data
    }

def define_temp_response_biochem():
    """Mocks define_temp_response_biochem.m"""
    # Simplified mock for the TDP structure
    return {'delHaV': 70000, 'delSV': 700, 'delHdV': 200000, 'delHaR': 46390, 'delSR': 100, 'delHdR': 200000, 'delHaKc': 79430, 'delHaKo': 36380, 'delHaT': 37830}

def assignvarnames():
    """Mocks assignvarnames.m"""
    # Simplified structure to match the variable names and expected outputs
    names = ['lite', 'calc_fluor', 'calc_planck', 'calc_xanthophyllabs', 'soilspectrum', 'Fluorescence_model', 'apply_T_corr', 'verify', 'saveCSV', 'mSCOPE', 'simulation', 'calc_directional', 'calc_vert_profiles', 'soil_heat_method', 'calc_rss_rbs', 'MoninObukhov', 'save_spectral', 'Simulation_Name', 'soil_file', 'optipar_file', 'atmos_file', 'Dataset_dir', 'meteo_ec_csv', 'vegetation_retrieved_csv', 'LIDF_file', 'verification_dir', 'mSCOPE_csv', 'nly', 't', 'year', 'Rin', 'Rli', 'p', 'Ta', 'ea', 'u', 'RH', 'VPD', 'tts', 'tto', 'psi', 'Cab', 'Cca', 'Cdm', 'Cw', 'Cs', 'Cant', 'N', 'SMC', 'BSMBrightness', 'BSMlat', 'BSMlon', 'LAI', 'hc', 'LIDFa', 'LIDFb', 'z', 'Ca', 'Vcmax25', 'BallBerrySlope', 'fqe', 'atmos_names']
    return [{'Name': name, 'Val': -999} for name in names]

def load_timeseries(V, F, xyt_in, path_input):
    """Mocks load_timeseries.m"""
    # Mocking a time series of 5 points
    telmax = 5
    t_mock = np.arange(telmax) * 3600 + datetime.now().timestamp()
    xyt = {'t': t_mock, 'year': np.zeros(telmax)}
    
    # Mock V update (only p, Ta, ea, Rin, Rli are time-dependent)
    for i, name in enumerate(['p', 'Ta', 'ea', 'Rin', 'Rli']):
        idx = next(i for i, v in enumerate(V) if v['Name'] == name)
        V[idx]['Val'] = np.ones(telmax) * V[idx]['Val']
        
    mly_ts = {}
    atmo_paths = [os.path.join(path_input, 'radiationdata', F[3].get('FileName', 'default_atmo_file.txt'))] * telmax
    return V, xyt, mly_ts, atmo_paths

def select_input(V, vi, canopy_in, options, constants, xyt_in=None, soil_in=None, leafbio_in=None):
    """Mocks select_input.m"""
    # This is a critical function that pulls variables from V based on index vi
    # and populates the model structures (soil, leafbio, canopy, meteo, angles, xyt).
    
    # Simple mocking for the structures:
    canopy = canopy_in.copy()
    
    # Mock angle and meteo selection based on index vi[tts/tto]
    angles = {'tts': 30, 'tto': 0, 'psi': 0}
    meteo = {'Ta': 298, 'ea': 20, 'p': 1013, 'u': 2, 'Rin': 800, 'Rli': 300, 'Ca': 400, 'z': 20}
    
    # Mock structure updates based on V and vi
    LAI_idx = next(i for i, v in enumerate(V) if v['Name'] == 'LAI')
    canopy['LAI'] = V[LAI_idx]['Val'][vi[LAI_idx]-1] if V[LAI_idx]['Val'].ndim > 0 else V[LAI_idx]['Val']
    
    # Mock xyt update
    if xyt_in is not None and 't' in xyt_in:
        t_idx = vi[next(i for i, v in enumerate(V) if v['Name'] == 't')] - 1
        xyt = {'t': xyt_in['t'][t_idx:t_idx+1], 'year': xyt_in['year'][t_idx:t_idx+1]}
    else:
        xyt = {'t': np.array([0]), 'year': np.array([0])}

    # Mock soil/leafbio structure selection
    soil = soil_in.copy() if soil_in is not None else {}
    leafbio = leafbio_in.copy() if leafbio_in is not None else {}

    # Update meteo time-series variables
    if options['simulation'] == 1:
        for name in ['p', 'Ta', 'ea', 'u', 'Rin', 'Rli']:
            idx = next(i for i, v in enumerate(V) if v['Name'] == name)
            k_val = vi[idx] - 1
            meteo[name] = V[idx]['Val'][k_val]

    # Mock Leafbio parameters
    for name in ['Cab', 'Cca', 'Cdm', 'Cw', 'Cs', 'Cant', 'N', 'Vcmax25', 'BallBerrySlope', 'fqe']:
        idx = next(i for i, v in enumerate(V) if v['Name'] == name)
        leafbio[name] = V[idx]['Val'][0] if V[idx]['Val'].ndim > 0 else V[idx]['Val']

    leafbio['rho_thermal'] = 0.01 # Placeholder for FLUSPECT
    leafbio['tau_thermal'] = 0.01 # Placeholder for FLUSPECT

    return soil, leafbio, canopy, meteo, angles, xyt

def load_atmo(atmfile, scopespec):
    """Mocks load_atmo.m"""
    # Mocking atmospheric data required by RTMo
    atmo = {}
    atmo['Esun_'] = np.ones(1800) * 1e-4 # Mocking Solar Irradiance
    atmo['Esky_'] = np.ones(1800) * 1e-5 # Mocking Sky Irradiance
    return atmo

def create_output_files_binary(parameter_file, F, path_of_code, path_input, spectral, options):
    """Mocks create_output_files_binary.m"""
    # Mocking output directory and file names
    Output_dir = os.path.join(path_of_code, 'output_mock')
    os.makedirs(Output_dir, exist_ok=True)
    f = {'t': 1, 'y': 2, 'R': 3} # Mock file handles/IDs
    fnames = ['output.bin']
    return Output_dir, f, fnames

def leafangles(LIDFa, LIDFb):
    """Mocks leafangles (ladgen).m"""
    # Simple placeholder for Leaf Inclination Distribution Function
    return np.ones((13, 36)) * (1 / (13 * 36))

def fluspect_mSCOPE(mly, spectral, leafbio, optipar, nl):
    """Mocks fluspect_mSCOPE.m"""
    # Mocking FLUSPECT RTM output
    nwl = 1800 # Mock number of wavelengths
    leafopt = {}
    leafopt['refl'] = np.ones((nl, nwl)) * 0.1
    leafopt['tran'] = np.ones((nl, nwl)) * 0.1
    return leafopt

def BSM(soil, optipar, soilemp):
    """Mocks BSM.m"""
    # Mocking BSM soil reflectance model
    return np.ones(1800) * 0.1

def RTMo(spectral, atmo, soil, leafopt, canopy, angles, constants, meteo, options):
    """Mocks RTMo.m (SAIL RTM)"""
    # Mocking RTMo outputs
    nl = canopy['nlayers']
    rad = {'Pnu_Cab': np.ones(nl) * 200, 'Pnh_Cab': np.ones(nl) * 100, 'Rnu_Cab': np.ones(nl) * 10, 'Rnh_Cab': np.ones(nl) * 5, 'Pnu_Car': np.zeros(nl), 'Pnh_Car': np.zeros(nl), 'Rnu_Car': np.zeros(nl), 'Rnh_Car': np.zeros(nl), 'Pnu': np.zeros(nl), 'Pnh': np.zeros(nl), 'Rnu_PAR': np.zeros(nl), 'Rnh_PAR': np.zeros(nl)}
    gap = {'Ps': np.ones(nl + 1) * 0.5}
    profiles = {}
    
    # Add other required rad fields for later steps
    rad['Lot'] = 1.0; rad['Lote'] = 1.0; rad['Lo_'] = np.zeros(1800); rad['Lot_'] = np.zeros(1800); rad['Eout_'] = np.zeros(1800); rad['Eoutte_'] = np.zeros(1800); rad['refl'] = np.zeros(1800)
    
    return rad, gap, profiles

def RTMf(constants, spectral, rad, soil, leafopt, canopy, gap, angles, bcu_eta, bch_eta):
    """Mocks RTMf.m (Fluorescence RTM)"""
    # Mocking RTMf outputs
    rad['LoF_'] = np.ones(spectral['wlF'].size) * 0.01
    return rad

def RTMz(constants, spectral, rad, soil, leafopt, canopy, gap, angles, bcu_Kn, bch_Kn):
    """Mocks RTMz.m (PRI RTM)"""
    # Mocking RTMz outputs
    return rad

def RTMt_planck(spectral, rad, soil, leafopt, canopy, gap, Tcu, Tch, Tsu, Tsh):
    """Mocks RTMt_planck.m"""
    return rad

def ephoton(wavelengths, constants):
    """Mocks ephoton: E = h * c / lambda"""
    h, c = 6.62607015e-34, 299792458
    # E [J/photon]
    E_photon = h * c / wavelengths 
    # E [J/umol] = E_photon * Avogadro's number
    return E_photon * constants['A'] * 1e-6

def Sint(spectrum, wavelength_axis):
    """Mocks Sint (Spectral Integration): Trapezoidal rule approximation"""
    return np.trapz(spectrum, wavelength_axis)

def calc_brdf(constants, options, directional, spectral, angles, atmo, soil, leafopt, canopy, meteo, thermal, bcu, bch):
    """Mocks calc_brdf.m"""
    return directional

def savebrdfoutput(options, directional, angles, spectral, Output_dir):
    """Mocks savebrdfoutput.m"""
    pass

def output_data_binary(f, k, xyt, rad, canopy, V, vi, vmax, options, fluxes, meteo, iter_data, resistance):
    """Mocks output_data_binary.m"""
    # Mocking a valid output process, returns column count
    return 100

def fill_output_with_nans(canopy, spectral):
    """Mocks fill_output_with_nans.m"""
    return canopy, {}, {}, {}, {}

def count_k(nvars, vi, vmax, step):
    """Mocks count_k for look-up table iteration"""
    # Simple increment for the first variable only
    vi[0] = vi[0] + step
    return vi

def bin_to_csv(fnames, V, vmax, n_col, telmax):
    """Mocks bin_to_csv.m"""
    pass

def output_verification_csv(Output_dir, verification_dir):
    """Mocks output_verification_csv.m"""
    pass

def input_mSCOPE(mSCOPE_csv_path):
    """Mocks input_mSCOPE.m"""
    return {'nly': 1, 'pLAI': 1.0, 'totLAI': 1.0, 'pCab': 50, 'pCca': 10, 'pCdm': 0.01, 'pCw': 0.01, 'pCs': 0.01, 'pN': 1.5}

# --- Main Execution Block ---

def scope_forward_simulation():
    """Main execution function translating SCOPE.m"""
    
    print('SCOPE forward simulation started')

    # Note: In Python, restoredefaultpath, addpath, clear all are not necessary
    # The current working directory is implicitly the starting path
    
    # 1. Define constants
    constants = define_constants()

    # 2. Paths
    path_input = 'input/'
    path_of_code = os.getcwd() # MATLAB's cd
    
    # --- 3. Simulation Options (Simplified Reading) ---
    
    # Mock reading the parameter filenames and N array for options
    # N is read from the first file, parameter file names from the second.
    # In a real setup, this would use a robust CSV reader.
    N_mock = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 2, 0, 0, 0] # Mock options values
    
    options = {}
    options['lite'] = N_mock[0]
    options['calc_fluor'] = N_mock[1]
    options['calc_planck'] = N_mock[2]
    options['calc_xanthophyllabs'] = N_mock[3]
    options['soilspectrum'] = N_mock[4]
    options['Fluorescence_model'] = N_mock[5]
    options['apply_T_corr'] = N_mock[6]
    options['verify'] = N_mock[7]
    options['saveCSV'] = N_mock[8]
    options['mSCOPE'] = N_mock[9]
    options['simulation'] = N_mock[10]
    options['calc_directional'] = N_mock[11]
    options['calc_vert_profiles'] = N_mock[12]
    options['soil_heat_method'] = N_mock[13]
    options['calc_rss_rbs'] = N_mock[14]
    options['MoninObukhov'] = N_mock[15]
    options['save_spectral'] = N_mock[16]
    
    if not (0 <= options['simulation'] <= 3):
        print('\n simulation option should be between 0 and 2 \r'); return
        
    options['Cca_function_of_Cab'] = 0
    
    integr = 'layers' if options['lite'] != 0 else 'angles_and_layers'

    # --- 3. File Names (F structure) ---
    f_names_list = ['Simulation_Name','soil_file','optipar_file','atmos_file', 'Dataset_dir', 'meteo_ec_csv', 'vegetation_retrieved_csv', 'LIDF_file', 'verification_dir', 'mSCOPE_csv', 'nly']
    cols_list = ['t', 'year', 'Rin','Rli', 'p','Ta','ea','u','RH', 'VPD', 'tts','tto', 'psi', 'Cab','Cca','Cdm','Cw','Cs','Cant','N', 'SMC','BSMBrightness', 'BSMlat', 'BSMlon', 'LAI', 'hc', 'LIDFa', 'LIDFb', 'z','Ca', 'Vcmax25', 'BallBerrySlope', 'fqe', 'atmos_names']
    fnc = f_names_list + cols_list
    
    # F = struct('FileID', fnc)
    F = [{'FileID': name, 'FileName': None} for name in fnc]
    
    # Mock file names assigned from the second parameter file
    mock_filenames = {
        'soil_file': 'soil_spectrum.csv', 
        'optipar_file': 'optipar.txt', 
        'atmos_file': 'modtran_default.txt', 
        'LIDF_file': 'LIDF_flat.csv'
    }
    for item in F:
        if item['FileID'] in mock_filenames:
            item['FileName'] = mock_filenames[item['FileID']]
    
    # --- 4. Input Data (V structure) ---
    
    # Mock reading input data from the third parameter file
    V = assignvarnames()
    
    # Mocking some input values (non-time series)
    input_values_mock = {
        'LAI': 2.5, 'hc': 2.0, 'z': 20.0, 'Ca': 400.0,
        'Cab': 40.0, 'Cca': 10.0, 'N': 1.4, 'Vcmax25': 60.0, 
        'BallBerrySlope': 9.0, 'Rin': 800.0, 'Rli': 300.0, 'p': 1013, 'Ta': 298.15
    }
    for item in V:
        if item['Name'] in input_values_mock:
            item['Val'] = np.array([input_values_mock[item['Name']]])

    # --- 6. Load Spectral Data ---
    
    # Mock load([path_input,'fluspect_parameters/', F(3).FileName])
    optipar = {'phi': np.ones(1800) * 0.01}
    
    # Mock soil spectrum load
    rsfile = np.ones((1800, 2)) * 0.1
    
    # --- 8. Define Canopy Structure ---
    canopy = {}
    canopy['nlincl'], canopy['nlazi'] = 13, 36
    canopy['litab'] = np.array([ 5, 15, 25, 35, 45, 55, 65, 75, 81, 83, 85, 87, 89])
    canopy['lazitab'] = np.arange(5, 360, 10)
    
    soilemp = {'SMC': 25, 'film': 0.015}
    
    LIDF_file = next(item for item in F if item['FileID'] == 'LIDF_file')['FileName']
    
    # --- 10. Define Spectral Regions ---
    spectral = define_bands()

    # --- 11. Load Time Series Data ---
    soil = {}
    xyt = {}
    if options['simulation'] == 1:
        vi = np.ones(len(V), dtype=int)
        soil, leafbio_mock, canopy_mock, meteo_mock, angles_mock, xyt_mock = select_input(V, vi, canopy, options, constants)
        V, xyt, mly_ts, atmo_paths = load_timeseries(V, F, xyt_mock, path_input)
        
        # Re-initialize structures with time-series dimensions
        soil, leafbio, canopy, meteo, angles, xyt = select_input(V, vi, canopy, options, constants)
    else:
        # Initial call for non-time series to populate structures
        vi = np.ones(len(V), dtype=int)
        soil, leafbio, canopy, meteo, angles, xyt = select_input(V, vi, canopy, options, constants)
        mly_ts = {}

    # --- 12. Preparations ---
    
    # Soil heat
    if options['simulation'] == 1 and options['soil_heat_method'] < 2:
        soil['Tsold'] = meteo['Ta'] * np.ones((12, 2))
        
    # Temperature sensitivity
    leafbio['TDP'] = define_temp_response_biochem()

    # Variables for loop
    vmax = np.array([len(v['Val']) for v in V])
    vmax[np.array([v['Name'] for v in V]) == 'LIDFa'] = 1 
    vmax[np.array([v['Name'] for v in V]) == 'LIDFb'] = 1 
    vi = np.ones(len(V), dtype=int)
    
    telmax = 0
    if options['simulation'] == 0: 
        telmax = np.max(vmax)
    elif options['simulation'] == 1: 
        telmax = xyt['t'].size
    elif options['simulation'] == 2:
        telmax = np.prod(vmax)
        
    xyt['t'], xyt['year'] = np.zeros(telmax), np.zeros(telmax)

    # Directional (BRDF)
    directional = {}
    if options['calc_directional']:
        anglesfile = np.array([[0, 0], [10, 0], [20, 0]]) # Mock angles
        directional['tto'] = anglesfile[:, 0]
        directional['psi'] = anglesfile[:, 1]
        directional['noa'] = len(directional['tto'])
    else:
        directional = None

    # Irradiance
    atmfile = os.path.join(path_input, 'radiationdata', next(item for item in F if item['FileID'] == 'atmos_file')['FileName'])
    atmo = load_atmo(atmfile, spectral['SCOPEspec'])

    # --- 13. Create output files ---
    Output_dir, f, fnames = create_output_files_binary(None, F, path_of_code, path_input, spectral, options)
    
    # --- 14. Run the models ---
    print('\n The calculations start now \r')
    
    start_time = datetime.now()
    
    for k in range(telmax):
        k_matlab = k + 1 # MATLAB index starts at 1
        
        # Update iteration index (vi)
        if options['simulation'] == 1: 
            vi[vmax > 1] = k_matlab
        if options['simulation'] == 0: 
            vi[vmax == telmax] = k_matlab

        # Select input for current time/run
        soil, leafbio, canopy, meteo, angles, xyt = select_input(V, vi, canopy, options, constants, xyt, soil, leafbio)
        
        # Canopy layers
        LAI = canopy['LAI']
        canopy['nlayers'] = math.ceil(10 * LAI) + ((meteo['Rin'] < 200) and options['MoninObukhov']) * 60
        canopy['nlayers'] = max(2, canopy['nlayers'])
        nl = canopy['nlayers']
        
        # Layer depths (x: from 0 to -1)
        x_col = np.arange(-1 / nl, -1 - 1e-9, -1 / nl)
        canopy['xl'] = np.insert(x_col, 0, 0) # [0; x]
        
        calculate = True
        
        if options['simulation'] != 1:
            print(f"simulation {k_matlab} of {telmax}")
        else:
            # Check for NaN in meteo inputs
            meteo_valid = not np.isnan(meteo['p'] * meteo['Ta'] * meteo['ea'] * meteo['u'] * meteo['Rin'] * meteo['Rli'])
            xyt_t_str = datetime.fromtimestamp(xyt['t'][0]).strftime('%Y-%m-%d %H:%M:%S')
            print(f"time = {xyt_t_str}: {k_matlab} / {telmax}")
            
            if not meteo_valid:
                canopy, fluxes, rad, resistance, iter_data = fill_output_with_nans(canopy, spectral)
                print('warning: run is invalid: there is NaN somewhere in meteo input')
                calculate = False
        
        if calculate:
            
            # LIDF (leafangles)
            if not LIDF_file:
                canopy['lidf'] = leafangles(canopy['LIDFa'], canopy['LIDFb'])
                
            # Leaf Radiative Transfer (FLUSPECT)
            leafbio['emis'] = 1 - leafbio['rho_thermal'] - leafbio['tau_thermal']
            leafbio['V2Z'] = 0
            
            # mSCOPE logic for vertical Cab/Cca profiles (simplified mock)
            mly = {}
            if options['simulation'] == 1 and mly_ts:
                # Use time series mly data (mocked to single values for simplicity)
                mly['nly'] = mly_ts['nly']
                mly['pLAI'] = mly_ts['pLAI'][k, :]
            else:
                mly['nly'] = 1
                mly['pLAI'] = LAI
                mly['totLAI'] = LAI
                mly['pCab'] = leafbio['Cab']
                mly['pCca'] = leafbio['Cca']
                mly['pCdm'] = leafbio['Cdm']
                mly['pCw'] = leafbio['Cw']
                mly['pCs'] = leafbio['Cs']
                mly['pN'] = leafbio['N']
            
            # Load atmo if new file
            atmo_paths_available = options['simulation'] == 1 and 'atmo_paths' in locals() and atmo_paths
            if atmo_paths_available and k_matlab > 1 and atmo_paths[k] != atmo_paths[k - 1]:
                atmo = load_atmo(atmo_paths[k], spectral['SCOPEspec'])
                
            leafopt = fluspect_mSCOPE(mly, spectral, leafbio, optipar, nl)
            # Assign thermal properties
            leafopt['refl'][:, spectral['IwlT']] = leafbio['rho_thermal']
            leafopt['tran'][:, spectral['IwlT']] = leafbio['tau_thermal']
            
            # Xanthophyll Cycle (RTMz requires this)
            if options['calc_xanthophyllabs']:
                leafbio['V2Z'] = 1
                leafoptZ = fluspect_mSCOPE(mly, spectral, leafbio, optipar, nl)
                leafopt['reflZ'] = leafopt['refl'].copy()
                leafopt['tranZ'] = leafopt['tran'].copy()
                # Update P bands with 'Z' properties
                leafopt['reflZ'][:, spectral['IwlP']] = leafoptZ['refl'][:, spectral['IwlP']]
                leafopt['tranZ'][:, spectral['IwlP']] = leafoptZ['tran'][:, spectral['IwlP']]

            # Soil Reflectance
            if options['soilspectrum'] == 0:
                soil['refl'] = rsfile[:, soil.get('spectrum', 0) + 1] # assuming spectrum index is 0 if not provided
            else:
                soil['refl'] = BSM(soil, optipar, soilemp)
            soil['refl'][spectral['IwlT']] = soil.get('rs_thermal', 0.99)

            # RTMo (incident radiation)
            rad, gap, profiles = RTMo(spectral, atmo, soil, leafopt, canopy, angles, constants, meteo, options)

            # Energy Balance (ebal)
            # Note: k_matlab is passed for time-step indexing in soil heat calculation
            iter_data, rad, thermal, soil, bcu, bch, fluxes, resistance, meteo = ebal(constants, options, rad, gap, meteo, soil, canopy, leafbio, k_matlab, xyt, integr)

            # Fluorescence RTM (RTMf)
            if options['calc_fluor']:
                rad = RTMf(constants, spectral, rad, soil, leafopt, canopy, gap, angles, bcu['eta'], bch['eta'])

            # PRI RTM (RTMz)
            if options['calc_xanthophyllabs']:
                # The corrected version of the RTMz call is:
                rad = RTMz(constants, spectral, rad, soil, leafopt, canopy, gap, angles, bcu['Kn'], bch['Kn'])
                
            # RTMt_sb (Final run to update thermal fields)
            rad = RTMt_sb(constants, rad, soil, leafbio, canopy, gap, thermal['Tcu'], thermal['Tch'], thermal['Tsu'], thermal['Tsh'], 1, spectral)
            
            # Planck RTM (calc_planck)
            if options['calc_planck']:
                rad = RTMt_planck(spectral, rad, soil, leafopt, canopy, gap, thermal['Tcu'], thermal['Tch'], thermal['Tsu'], thermal['Tsh'])

            # --- Data Products (LST, A, Ja, NPQ, SIF) ---
            
            Ps = gap['Ps'][:nl]
            Ph = (1 - Ps)

            # Canopy LAI fractions
            canopy['LAIsunlit'] = LAI * np.mean(Ps)
            canopy['LAIshaded'] = LAI - canopy['LAIsunlit']

            # Aggregated PAR/Radiance fluxes (Pntot, Rntot) - requires meanleaf
            canopy['Pnsun_Cab'] = LAI * meanleaf(canopy, rad['Pnu_Cab'], integr, Ps)
            canopy['Pnsha_Cab'] = LAI * meanleaf(canopy, rad['Pnh_Cab'], 'layers', Ph)
            canopy['Pntot_Cab'] = canopy['Pnsun_Cab'] + canopy['Pnsha_Cab']

            # LST (black-body surface assumption)
            canopy['LST'] = (np.pi * (rad['Lot'] + rad['Lote']) / (constants['sigmaSB'] * rad['canopyemis']))**0.25
            canopy['emis'] = rad['canopyemis']

            # Photosynthesis
            canopy['A'] = LAI * (meanleaf(canopy, bch['A'], 'layers', Ph) + meanleaf(canopy, bcu['A'], integr, Ps))
            canopy['GPP'] = LAI * (meanleaf(canopy, bch['Ag'], 'layers', Ph) + meanleaf(canopy, bcu['Ag'], integr, Ps))
            canopy['Ja'] = LAI * (meanleaf(canopy, bch['Ja'], 'layers', Ph) + meanleaf(canopy, bcu['Ja'], integr, Ps))
            
            # NPQ (Energy/Photons)
            canopy['PNPQ'] = LAI * (meanleaf(canopy, rad['Pnh_Cab'] * bch['Phi_N'], 'layers', Ph) + meanleaf(canopy, rad['Pnu_Cab'] * bcu['Phi_N'], integr, Ps))

            # Fluorescence Re-absorption corrected (SIF-reabsorption correction)
            aPAR_Cab_eta = LAI * (meanleaf(canopy, bch['eta'] * rad['Pnh_Cab'], 'layers', Ph) + meanleaf(canopy, bcu['eta'] * rad['Pnu_Cab'], integr, Ps)) 

            if options['calc_fluor']:
                # Constants for fluorescence
                ep = constants['A'] * ephoton(spectral['wlF'] * 1e-9, constants)
                rad['PoutFrc'] = leafbio['fqe'] * aPAR_Cab_eta
                # EoutFrc_ = 1E-3*ep.*(rad.PoutFrc*optipar.phi(spectral.IwlF))
                # Note: optipar.phi(spectral.IwlF) is a scaling factor for the spectrum
                IwlF_size = spectral['wlF'].size
                rad['EoutFrc_'] = 1e-3 * ep * (rad['PoutFrc'] * optipar['phi'][:IwlF_size])
                
                rad['EoutFrc'] = 1e-3 * Sint(rad['EoutFrc_'], spectral['wlF'])
                
                # Fluorescence yield
                sigmaF = np.pi * rad['LoF_'] / rad['EoutFrc_']
                rad['sigmaF'] = np.interp(spectral['wlF'], spectral['wlF'], sigmaF) # Simplified interp
                
                canopy['fqe'] = rad['PoutFrc'] / canopy['Pntot_Cab']
            else:
                canopy['fqe'] = np.nan

            # Final composite radiance/reflectance
            rad['Lotot_'] = rad['Lo_'] + rad['Lot_']
            rad['Eout_'] = rad['Eout_'] + rad['Eoutte_']
            
            if options['calc_fluor']:
                rad['Lototf_'] = rad['Lotot_'].copy()
                IwlF_indices = np.where((spectral['wl'] >= spectral['wlF'][0]) & (spectral['wl'] <= spectral['wlF'][-1]))[0]
                rad['Lototf_'][IwlF_indices] = rad['Lototf_'][IwlF_indices] + rad['LoF_']
                rad['reflapp'] = rad['refl'].copy()
                rad['reflapp'][IwlF_indices] = np.pi * rad['Lototf_'][IwlF_indices] / (atmo['Esun_'][IwlF_indices] + atmo['Esky_'][IwlF_indices])
            
            # Directional output
            if options['calc_directional']:
                directional = calc_brdf(constants, options, directional, spectral, angles, atmo, soil, leafopt, canopy, meteo, thermal, bcu, bch)
                savebrdfoutput(options, directional, angles, spectral, Output_dir)

            # Integrated Lo
            IwlP_indices = spectral['IwlP']
            rad['Lo'] = 0.001 * Sint(rad['Lo_'][IwlP_indices], spectral['wl'][IwlP_indices])
            
            # --- Write output ---
            n_col = output_data_binary(f, k_matlab, xyt, rad, canopy, V, vi, vmax, options, fluxes, meteo, iter_data, resistance)

            # --- Update input (Look-up table) ---
            if options['simulation'] == 2 and telmax > 1:
                vi = count_k(len(V), vi, vmax, 1)

    # Final steps
    end_time = datetime.now()
    print(f"Time elapsed: {end_time - start_time}")
    
    # fclose('all') is implicit in Python
    
    if options['saveCSV']:
        bin_to_csv(fnames, V, vmax, n_col, telmax)
        
    if options['verify']:
        output_verification_csv(Output_dir, next(item for item in F if item['FileID'] == 'verification_dir')['FileName'])

# Run the main function
if __name__ == '__main__':
    scope_forward_simulation()