import numpy as np
import struct
import warnings

# Assume helper functions are available
from .timestamp_utils import juliandate, timestamp2datetime # Mock/Placeholder for these

# --- Placeholder/Mock Functions ---

# Placeholder for Julian date calculation (juliandate)
def juliandate(dt):
    """MOCK: Converts datetime/timestamp to Julian date number."""
    # This is highly dependent on the type of dt. Assume conversion is correct.
    return dt.toordinal()

# --- Main Function ---

def output_data_binary(f, k, xyt, rad, canopy, V, vi, vmax, options, fluxes, meteo, iter_out, resistance):
    """
    Translation of src/IO/output_data_binary.m
    Writes output data for a single timestep 'k' to binary files.
    """
    
    # --- Time Conversion (for output) ---
    # MATLAB: Convert internal datetime to DOY for output if needed
    if isinstance(xyt['t'][k], pd.Timestamp):
        # We need the full series for the start/end DOY conversion
        t_full_dt = pd.to_datetime(xyt['t'])
        
        # MATLAB: get_doy = @(x) juliandate(x) - juliandate(datetime(year(x), 1, 0));
        # Simplified DOY calculation for Pandas Timestamp (1-based)
        doy_series = t_full_dt.dt.dayofyear
        
        # Update V with start/end DOY (V indices 46 and 47 are 45 and 46)
        if V[45]['Name'] == 'startDOY':
            V[45]['Val'] = doy_series.min()
        if V[46]['Name'] == 'endDOY':
            V[46]['Val'] = doy_series.max()

        # Current time step DOY
        current_doy = doy_series.iloc[k]
    else:
        # Assumed to be DOY or numerical already
        current_doy = xyt['t'][k]
        
    # --- Prepare for writing (all data is written as 'double', 8 bytes) ---
    def write_double(file_handle, data):
        # Ensure data is flattened and a NumPy array for consistent writing
        data_flat = np.array(data).flatten()
        # MATLAB's 'double' is 8-byte float. Python's struct uses 'd'.
        packed_data = struct.pack(f'<{len(data_flat)}d', *data_flat)
        file_handle.write(packed_data)
        return len(data_flat)

    n_col = {} # Stores column counts

    # %% Vegetation products
    apar_out = [
        k + 1, xyt['year'][k], current_doy, rad['PAR'], rad['EPAR'], canopy['LAIsunlit'], canopy['LAIshaded'],
        canopy['Pntot'], canopy['Pnsun'], canopy['Pnsha'],
        canopy['Pntot_Cab'], canopy['Pnsun_Cab'], canopy['Pnsha_Cab'],
        canopy['Pntot_Car'], canopy['Pnsun_Car'], canopy['Pnsha_Car'],
        canopy['Rntot_PAR'], canopy['Rnsun_PAR'], canopy['Rnsha_PAR'],
        canopy['Rntot_Cab'], canopy['Rnsun_Cab'], canopy['Rnsha_Cab'],
        canopy['Rntot_Car'], canopy['Rnsun_Car'], canopy['Rnsha_Car']
    ]
    n_col['apar'] = write_double(f['apar_file'], apar_out)

    veg_out = [
        k + 1, xyt['year'][k], current_doy, canopy['A'], canopy['Ja'], canopy['ENPQ'], 
        canopy['PNPQ'], canopy['fqe'], canopy['LST'], canopy['emis'], canopy['GPP']
    ]
    n_col['veg'] = write_double(f['veg_file'], veg_out)

    # %% Fluxes product
    # MATLAB: cell2mat(struct2cell(fluxes))'
    # Need to preserve the order of fields in fluxes
    fluxes_ordered = [
        fluxes.get('Rnctot', np.nan), fluxes.get('lEctot', np.nan), fluxes.get('Hctot', np.nan), fluxes.get('Actot', np.nan), fluxes.get('Tcave', np.nan),
        fluxes.get('Rnstot', np.nan), fluxes.get('lEstot', np.nan), fluxes.get('Hstot', np.nan), fluxes.get('Gtot', np.nan), fluxes.get('Tsave', np.nan),
        fluxes.get('Rntot', np.nan), fluxes.get('lEtot', np.nan), fluxes.get('Htot', np.nan)
    ]
    flu_out = [k + 1, iter_out['counter'], xyt['year'][k], current_doy] + fluxes_ordered
    n_col['flu'] = write_double(f['flu_file'], flu_out)

    # %% Radiation
    rad_out = [
        k + 1, xyt['year'][k], current_doy, meteo['Rin'], meteo['Rli'], rad['Eouto'], 
        rad['Eoutt'] + rad['Eoutte'], rad['Lo'], rad['Lot'], rad['Lote']
    ]
    n_col['rad'] = write_double(f['rad_file'], rad_out)

    # %% Fluorescence spectral and scalar outputs
    if options['calc_fluor']:
        fluor_out = [
            rad['F685'], rad['wl685'], rad['F740'], rad['wl740'], rad['F684'], rad['F761'], 
            rad['LoutF'], rad['EoutF'], rad['EoutFrc']
        ]
        n_col['fluor'] = write_double(f['fluor_file'], fluor_out)

        n_col['fluor_spectrum'] = write_double(f['fluor_spectrum_file'], rad['LoF_'])
        n_col['sigmaF'] = write_double(f['sigmaF_file'], rad['sigmaF'])
        n_col['fRC'] = write_double(f['fRC_file'], rad['EoutFrc_'])
        n_col['fRCL'] = write_double(f['fRCL_file'], rad['Femleaves_'])
        n_col['fhemis'] = write_double(f['fhemis_file'], rad['EoutF_'])
        n_col['Lo2'] = write_double(f['Lo2_file'], rad['Lototf_'])
        n_col['rapp'] = write_double(f['rapp_file'], rad['reflapp'])

    # %% reflectance and Radiance
    if 'r_file' in f:
        n_col['r'] = write_double(f['r_file'], rad['refl'])
        n_col['rsd'] = write_double(f['rsd_file'], rad['rsd'])
        n_col['rdd'] = write_double(f['rdd_file'], rad['rdd'])
        n_col['rso'] = write_double(f['rso_file'], rad['rso'])
        n_col['rdo'] = write_double(f['rdo_file'], rad['rdo'])

        n_col['Eout'] = write_double(f['Eout_file'], rad['Eout_'])
        n_col['Lo'] = write_double(f['Lo_file'], rad['Lotot_'])
        n_col['Esun'] = write_double(f['Esun_file'], rad['Esun_'])
        n_col['Esky'] = write_double(f['Esky_file'], rad['Esky_'])

    # %% Resistances
    resist_out = [resistance['raa'], resistance['raws'], resistance['rss'], resistance['ustar']]
    n_col['resist'] = write_double(f['resist_file'], resist_out)

    # %% pars (short parameter file)
    k2 = np.where(vmax > 1)[0]
    V_short = np.empty(len(k2) + 1)
    V_short[0] = len(k2) # number of parameters
    for i, idx in enumerate(k2):
        # V(idx).Val(vi(idx)) -> Val array, indexed by vi at current k
        V_short[i + 1] = V[idx]['Val'][vi[idx][k]]

    n_col['pars'] = write_double(f['pars_file'], V_short)

    return n_col