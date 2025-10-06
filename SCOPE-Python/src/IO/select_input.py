import numpy as np
import warnings

# Assume external helper functions are defined elsewhere
# from .soil_inertia import Soil_Inertia1, Soil_Inertia0
# from .calc_rssrbs import calc_rssrbs
# from .zo_and_d import zo_and_d

# --- Placeholder/Mock Functions ---

# Mock the derived input functions
def Soil_Inertia1(SMC): return np.mean(SMC) * 10 
def Soil_Inertia0(cs, rhos, lambdas): return 5.0
def calc_rssrbs(SMC, LAI, rbs): return 100.0, 50.0
def zo_and_d(soil, canopy, constants): return canopy['zo'], canopy['d']

# --- Main Function ---

def select_input(V, vi, canopy, options, constants, xyt, soil, leafbio):
    """
    Translation of src/IO/select_input.m
    Populates model structures with parameter values from the V list for the current timestep.
    """
    
    # Helper to get the value for the current timestep 'k'
    def get_val(v_idx, k):
        # v_idx is 1-based index in MATLAB, so idx-1 in Python
        # V[idx-1]['Val'] is the full time series; [k] gets the value at k
        return V[v_idx - 1]['Val'][k]

    k = vi[0] # The current time step index 'k' in the original MATLAB array 'vi' (which contained 1s or k)
    
    # The MATLAB implementation passes vi as a full array (length of V) where each element is the index
    # (usually k) to select from the time series, or 1 if it's a scalar.
    # We will assume vi is a list/array of indices (len(V)) and k is the index used to access the
    # current time step value for each V['Val'] array.
    
    # We iterate through V and use the value V[i]['Val'][vi[i]]
    # Since V is 0-indexed in Python, V(n).Val(vi(n)) becomes V[n-1]['Val'][vi[n-1]]
    
    # We must assume vi is an array of indices corresponding to the full V list.
    
    # --- Soil Parameters ---
    soil['spectrum'] = get_val(60, vi[59])
    soil['rss'] = get_val(64, vi[63])
    soil['rs_thermal'] = get_val(65, vi[64])
    soil['cs'] = get_val(66, vi[65])
    soil['rhos'] = get_val(67, vi[66])
    soil['CSSOIL'] = get_val(43, vi[42])
    soil['lambdas'] = get_val(68, vi[67])
    soil['rbs'] = get_val(44, vi[43])
    soil['SMC'] = get_val(54, vi[53])
    soil['BSMBrightness'] = get_val(61, vi[60])
    soil['BSMlat'] = get_val(62, vi[61])
    soil['BSMlon'] = get_val(63, vi[62])

    # --- Leaf Biochemical/Optical Parameters ---
    leafbio['Cab'] = get_val(1, vi[0])
    leafbio['Cca'] = get_val(2, vi[1])
    if options['Cca_function_of_Cab']:
        leafbio['Cca'] = 0.25 * leafbio['Cab']
    leafbio['Cdm'] = get_val(3, vi[2])
    leafbio['Cw'] = get_val(4, vi[3])
    leafbio['Cs'] = get_val(5, vi[4])
    leafbio['N'] = get_val(6, vi[5])
    leafbio['Cant'] = get_val(7, vi[6])
    leafbio['Cp'] = get_val(8, vi[7])
    leafbio['Cbc'] = get_val(9, vi[8])
    leafbio['rho_thermal'] = get_val(10, vi[9])
    leafbio['tau_thermal'] = get_val(11, vi[10])

    # --- Leaf Physiological Parameters ---
    leafbio['Vcmax25'] = get_val(12, vi[11])
    leafbio['BallBerrySlope'] = get_val(13, vi[12])
    leafbio['BallBerry0'] = get_val(14, vi[13])
    leafbio['Type'] = 'C4' if get_val(15, vi[14]) else 'C3' # MATLAB: Type is 1 for C4, 0 for C3 in input
    canopy['kV'] = get_val(16, vi[15])
    leafbio['RdPerVcmax25'] = get_val(17, vi[16])
    fqe_val = get_val(18, vi[17])
    leafbio['Kn0'] = get_val(19, vi[18])
    leafbio['Knalpha'] = get_val(20, vi[19])
    leafbio['Knbeta'] = get_val(21, vi[20])

    leafbio['Tyear'] = get_val(55, vi[54])
    leafbio['beta'] = get_val(56, vi[55])
    leafbio['kNPQs'] = get_val(57, vi[56])
    leafbio['qLs'] = get_val(58, vi[57])
    leafbio['stressfactor'] = get_val(59, vi[58])

    # --- Canopy Parameters ---
    canopy['LAI'] = max(1E-9, get_val(22, vi[21]))
    canopy['hc'] = get_val(23, vi[22])
    canopy['zo'] = get_val(24, vi[23])
    canopy['d'] = get_val(25, vi[24])
    canopy['LIDFa'] = get_val(26, vi[25])
    canopy['LIDFb'] = get_val(27, vi[26]) # Note: MATLAB uses index 26 twice here for LIDFb
    canopy['leafwidth'] = get_val(28, vi[27])
    canopy['rb'] = get_val(38, vi[37])
    canopy['Cd'] = get_val(39, vi[38])
    canopy['CR'] = get_val(40, vi[39])
    canopy['CD1'] = get_val(41, vi[40])
    canopy['Psicor'] = get_val(42, vi[41])
    canopy['rwc'] = get_val(45, vi[44])

    # --- Meteorological Inputs ---
    meteo['z'] = get_val(29, vi[28])
    meteo['Rin'] = get_val(30, vi[29])
    meteo['Ta'] = get_val(31, vi[30])
    meteo['Rli'] = get_val(32, vi[31])
    meteo['p'] = get_val(33, vi[32])
    meteo['ea'] = get_val(34, vi[33])
    meteo['u'] = get_val(35, vi[34])
    meteo['Ca'] = get_val(36, vi[35])
    meteo['Oa'] = get_val(37, vi[36])

    # --- Site/Angle/Time Parameters (xyt and angles) ---
    xyt['startDOY'] = get_val(46, vi[45])
    xyt['endDOY'] = get_val(47, vi[46])
    xyt['LAT'] = get_val(48, vi[47])
    xyt['LON'] = get_val(49, vi[48])
    xyt['timezn'] = get_val(50, vi[49])

    angles['tts'] = get_val(51, vi[50])
    angles['tto'] = get_val(52, vi[51])
    angles['psi'] = get_val(53, vi[52])

    # Final check on SMC units
    if soil['SMC'] > 1:
        soil['SMC'] /= 100

    # --- Derived Input Calculations ---
    if options['soil_heat_method'] == 1:
        soil['GAM'] = Soil_Inertia1(soil['SMC'])
    else:
        soil['GAM'] = Soil_Inertia0(soil['cs'], soil['rhos'], soil['lambdas'])

    if options['calc_rss_rbs']:
        soil['rss'], soil['rbs'] = calc_rssrbs(soil['SMC'], canopy['LAI'], soil['rbs'])

    # leafbio.Type is already set above via V(15).Val(vi(15))
    canopy['hot'] = canopy['leafwidth'] / canopy['hc']
    canopy['zo'], canopy['d'] = zo_and_d(soil, canopy, constants)
    leafbio['fqe'] = fqe_val # Assign the fqe value for consistency

    return soil, leafbio, canopy, meteo, angles, xyt