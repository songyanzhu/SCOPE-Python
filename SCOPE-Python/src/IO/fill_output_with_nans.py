import numpy as np

def fill_output_with_nans(canopy, spectral):
    """
    Translation of src/IO/fill_output_with_nans.m
    Initializes output structures with NaN values for the current simulation step.
    """

    # canopy is passed as an existing dictionary and is updated in place (implicit in MATLAB)
    
    # %% canopy
    canopy_names = [
        "LAIsunlit", "LAIshaded", "Pntot", "Pnsun", "Pnsha",
        "Pntot_Cab", "Pnsun_Cab", "Pnsha_Cab", "Pntot_Car", "Pnsun_Car", "Pnsha_Car",
        "Rntot_PAR", "Rnsun_PAR", "Rnsha_PAR", "Rntot_Cab", "Rnsun_Cab", "Rnsha_Cab",
        "Rntot_Car", "Rnsun_Car", "Rnsha_Car", "A", "Ja", "ENPQ", "PNPQ", "fqe", "LST", "emis"
    ]
    for name in canopy_names:
        canopy[name] = np.nan

    # %% fluxes
    fluxes_names = [
        'Rnctot','lEctot','Hctot','Actot','Tcave', 
        'Rnstot','lEstot','Hstot','Gtot','Tsave',
        'Rntot','lEtot','Htot'
    ]
    fluxes = {}
    for name in fluxes_names:
        fluxes[name] = np.nan

    # %% rad scalars
    rad_names_scalar = [
        "PAR", "EPAR", "Eouto","Eoutt","Eoutte", "Lo","Lot","Lote",
        "F685","wl685","F740","wl740","F684","F761", "LoutF","EoutF","EoutFrc"
    ]
    rad = {}
    for name in rad_names_scalar:
        rad[name] = np.nan

    # %% rad spectral (size of wlS)
    wlS_size = spectral['wlS'].shape
    rad_names_spectral = [
        "reflapp", "refl", "rsd", "rdd", "rso", "rdo", 
        "Eout_", "Lotot_", "Lototf_", "Esun_", "Esky_"
    ]
    for name in rad_names_spectral:
        rad[name] = np.full(wlS_size, np.nan)

    # %% sif spectral (size of wlF)
    wlF_size = spectral['wlF'].shape
    sif_names = [
        "LoF_", "sigmaF", "EoutFrc_", "Femleaves_", "EoutF_"
    ]
    for name in sif_names:
        rad[name] = np.full(wlF_size, np.nan)

    # %% resistances
    resistances_names = ['raa','raws','rss','ustar']
    resistance = {}
    for name in resistances_names:
        resistance[name] = np.nan

    # %% iter
    iter_out = {'counter': np.nan}
    
    # MATLAB returns the updated/new structures
    return canopy, fluxes, rad, resistance, iter_out