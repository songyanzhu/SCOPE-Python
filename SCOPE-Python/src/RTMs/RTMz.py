import numpy as np
# Assuming helper functions _ephoton, _e2phot, _sint, etc. are available in this module's scope
# as they were either defined in BSM.py or RTMf.py

# Placeholder for _ephoton, _e2phot, as they are not defined in this file.
def _ephoton(lambda_val, constants):
    return constants['h'] * constants['c'] / lambda_val

def _e2phot(lambda_val, E, constants):
    e = _ephoton(lambda_val, constants)
    return E / e / constants['A']

# Placeholder for Kn2Cx
def _kn2cx(Kn):
    """empirical fit by N Vilfan (Vilfan et al, 2018, 2019)"""
    return 0.3187 * Kn

# --- Main Function ---

def RTMz(constants, spectral, rad, soil, leafopt, canopy, gap, angles, Knu, Knh):
    """
    Translation of src/RTMs/RTMz.m
    Calculates the modification of TOC outgoing radiance due to the conversion 
    of Violaxanthin into Zeaxanthin (NPQ component).
    """

    # --- 0.1 initialisations ---
    
    wlS = spectral['wlS']
    wlZ = spectral['wlZ']
    
    # Intersect for indices
    iwlfi = np.isin(wlS, wlZ)
    iwlfi_idx = np.where(iwlfi)[0]
    nwlZ = len(wlZ)
    
    nl = canopy['nlayers']
    LAI, iLAI = canopy['LAI'], canopy['LAI'] / nl
    litab, laziab, lidf = canopy['litab'], canopy['lazitab'], canopy['lidf']
    nlazi, nlinc = canopy['nlazi'], canopy['nlinc']
    nlori = nlinc * nlazi
    layers = np.arange(1, nl + 1)

    Ps, Po, Pso = gap['Ps'], gap['Po'], gap['Pso']
    Qs = Ps[:-1] # Ps(1:end-1)

    # Difference in optical properties due to V2Z conversion
    RZ = (leafopt['reflZ'][:, iwlfi_idx] - leafopt['refl'][:, iwlfi_idx]).T # [nwlZ, nl]
    TZ = (leafopt['tranZ'][:, iwlfi_idx] - leafopt['tran'][:, iwlfi_idx]).T # [nwlZ, nl]

    # RTM coefficients at Z wavelengths
    Esunf_ = rad['Esun_'][iwlfi_idx] # [nwlZ]
    
    # Eminf_/Epluf_: Diffuse radiation components, [nwlZ, nl+1, 2 (sun/sky)]
    Eminf_ = np.zeros((nwlZ, nl + 1, 2))
    Epluf_ = np.zeros((nwlZ, nl + 1, 2))
    
    # Eminf_(:,:,1) = rad.Emins_(:,iwlfi)'
    Eminf_[:, :, 0] = rad['Emins_'][1:nl+2, iwlfi_idx].T
    Eminf_[:, :, 1] = rad['Emind_'][1:nl+2, iwlfi_idx].T
    Epluf_[:, :, 0] = rad['Eplus_'][:nl+1, iwlfi_idx].T
    Epluf_[:, :, 1] = rad['Eplud_'][:nl+1, iwlfi_idx].T

    # RTM coefficients at Z wavelengths [nl, nwlZ]
    Xdd = rad['Xdd'][:, iwlfi_idx]
    rho_dd = rad['rho_dd'][:, iwlfi_idx]
    R_dd = rad['R_dd'][:, iwlfi_idx]
    tau_dd = rad['tau_dd'][:, iwlfi_idx]
    vb = rad['vb'][:, iwlfi_idx]
    vf = rad['vf'][:, iwlfi_idx]

    # --- 0.2 geometric quantities (from RTMf) ---
    
    deg2rad = constants['deg2rad']
    tto, tts, psi = angles['tto'], angles['tts'], angles['psi']
    rs = soil['refl'][iwlfi_idx]
    
    cos_tto, sin_tto = np.cos(tto * deg2rad), np.sin(tto * deg2rad)
    cos_tts, sin_tts = np.cos(tts * deg2rad), np.sin(tts * deg2rad)
    cos_ttli, sin_ttli = np.cos(litab * deg2rad), np.sin(litab * deg2rad)
    cos_phils = np.cos(lazitab * deg2rad)
    cos_philo = np.cos((lazitab - psi) * deg2rad)

    # --- 0.3 geometric factors for all leaf angle/azumith classes (from RTMf) ---
    
    cos_tts_exp = np.outer(cos_ttli, np.ones(nlazi)) * cos_tts
    sin_tts_exp = np.outer(sin_ttli, np.ones(nlazi)) * sin_tts
    
    cds = cos_tts_exp + sin_tts_exp * cos_phils 
    cdo = np.outer(cos_ttli, np.ones(nlazi)) * cos_tto + np.outer(sin_ttli, np.ones(nlazi)) * sin_tto * cos_philo
    
    fs = cds / cos_tts
    absfs = np.abs(fs)
    fo = cdo / cos_tto
    absfo = np.abs(fo)
    fsfo = fs * fo
    absfsfo = np.abs(fsfo)
    
    cos_ttli_exp = np.outer(cos_ttli, np.ones(nlazi))
    foctl = fo * cos_ttli_exp
    fsctl = fs * cos_ttli_exp
    ctl2 = np.outer(cos_ttli**2, np.ones(nlazi))
    
    # reshape all the variables to [nlori, 1]
    absfs, absfo = absfs.reshape(-1, 1), absfo.reshape(-1, 1)
    fsfo, absfsfo = fsfo.reshape(-1, 1), absfsfo.reshape(-1, 1)
    foctl, fsctl = foctl.reshape(-1, 1), fsctl.reshape(-1, 1)
    ctl2 = ctl2.reshape(-1, 1)

    # --- 1.0 calculation of 'flux' in observation direction ---
    
    # Fmin_/Fplu_ is [nl+1, nwlZ, 2]
    Fmin_z, Fplu_z = np.zeros((nl + 1, nwlZ, 2)), np.zeros((nl + 1, nwlZ, 2))
    laz = 1 / nlazi
    
    etah = _kn2cx(Knh)

    if Knu.ndim == 1:
        etau = _kn2cx(Knu)
        etau = np.repeat(etau, nlinc * nlazi).reshape(nl, nlinc, nlazi)
        etau = np.transpose(etau, (1, 2, 0)) # [nlinc, nlazi, nl]
        etau_nlori_nl = etau.reshape(nlori, nl) # [nlori, nl]
    else:
        # Assuming Knu is [nl, nlinc, nlazi]
        etau = _kn2cx(Knu)
        etau_nlori_nl = np.transpose(etau, (1, 2, 0)).reshape(nlori, nl) 
    
    # etau_lidf [nlori, nl]
    lidf_laz = np.outer(lidf, np.ones(nlazi)) * laz
    lidf_laz_nlori = lidf_laz.reshape(nlori, 1)
    
    etau_lidf = etau_nlori_nl * lidf_laz_nlori
    etah_lidf = np.tile(etah[np.newaxis, :], (nlori, 1)) * lidf_laz_nlori

    # Loop for sun (k=1) and sky (k=2) incidence
    for k in range(2):
        U, Y = np.zeros((nl + 1, nwlZ)), np.zeros((nl, nwlZ))

        # Optical properties difference matrices (Mplu/Mmin are RZ/TZ here)
        MpluEsun = RZ * Esunf_[:, np.newaxis] * (k == 0) # k=0 is sun
        MminEsun = TZ * Esunf_[:, np.newaxis] * (k == 0)
        MpluEmin, MpluEplu = RZ * Eminf_[:, :nl, k], RZ * Epluf_[:, :nl, k]
        MminEmin, MminEplu = TZ * Eminf_[:, :nl, k], TZ * Epluf_[:, :nl, k]
        
        # All bsxfun operations follow RTMf structure (vectorized over nwlZ)
        
        # Emission factors (PiL)
        wfEs = MpluEsun * np.sum(etau_lidf * absfsfo, axis=0) + MminEsun * np.sum(etau_lidf * fsfo, axis=0)
        piLs = wfEs + MpluEplu * np.sum(etau_lidf * absfo, axis=0) - MminEplu * np.sum(etau_lidf * foctl, axis=0) + \
                      MpluEmin * np.sum(etau_lidf * absfo, axis=0) + MminEmin * np.sum(etau_lidf * foctl, axis=0)
        piLd = MpluEmin * np.sum(etah_lidf * absfo, axis=0) + MminEmin * np.sum(etah_lidf * foctl, axis=0) + \
               MpluEplu * np.sum(etah_lidf * absfo, axis=0) - MminEplu * np.sum(etah_lidf * foctl, axis=0)
        
        # Scattering factors (Fs/Fd)
        Fsmin = MpluEsun * np.sum(etau_lidf * absfs, axis=0) - MminEsun * np.sum(etau_lidf * fsctl, axis=0) + \
                MpluEmin * np.sum(etau_lidf, axis=0) - MminEmin * np.sum(etau_lidf * ctl2, axis=0) + \
                MpluEplu * np.sum(etau_lidf, axis=0) + MminEplu * np.sum(etau_lidf * ctl2, axis=0)
        
        Fsplu = MpluEsun * np.sum(etau_lidf * absfs, axis=0) + MminEsun * np.sum(etau_lidf * fsctl, axis=0) + \
                MpluEmin * np.sum(etau_lidf, axis=0) + MminEmin * np.sum(etau_lidf * ctl2, axis=0) + \
                MpluEplu * np.sum(etau_lidf, axis=0) - MminEplu * np.sum(etau_lidf * ctl2, axis=0)
        
        Fdmin = MpluEmin * np.sum(etah_lidf, axis=0) - MminEmin * np.sum(etah_lidf * ctl2, axis=0) + \
                MpluEplu * np.sum(etah_lidf, axis=0) + MminEplu * np.sum(etah_lidf * ctl2, axis=0)
        
        Fdplu = MpluEmin * np.sum(etah_lidf, axis=0) + MminEmin * np.sum(etah_lidf * ctl2, axis=0) + \
                MpluEplu * np.sum(etah_lidf, axis=0) - MminEplu * np.sum(etah_lidf * ctl2, axis=0)
        
        # Emitted RTM source terms
        Qs_exp = Qs[np.newaxis, :]
        Femmin = iLAI * (Fsmin * Qs_exp + Fdmin * (1 - Qs_exp))
        Femplu = iLAI * (Fsplu * Qs_exp + Fdplu * (1 - Qs_exp))
        
        # RTM (Bottom to Top)
        for j in range(nl - 1, -1, -1):
            rho_dd_j, R_dd_j1 = rho_dd[j], R_dd[j+1]
            Y[j] = (rho_dd_j * U[j+1] + Femmin[:, j]) / (1 - rho_dd_j * R_dd_j1)
            tau_dd_j = tau_dd[j]
            U[j] = tau_dd_j * (R_dd_j1 * Y[j] + U[j+1]) + Femplu[:, j]

        # RTM (Top to Bottom)
        for j in range(nl):
            Fmin_z[j+1, :, k] = Xdd[j] * Fmin_z[j, :, k] + Y[j]
            Fplu_z[j, :, k] = R_dd[j] * Fmin_z[j, :, k] + U[j]

        # Directional Radiance
        piLo1 = iLAI * piLs @ Pso[:nl]
        piLo2 = iLAI * piLd @ (Po[:nl] - Pso[:nl])
        
        Fmin_layers = Fmin_z[layers, :, k]
        Fplu_layers = Fplu_z[layers, :, k]
        Po_layers = Po[:nl]
        
        piLo3 = iLAI * (vb[:nl] * Fmin_layers + vf[:nl] * Fplu_layers).T @ Po_layers
        piLo4 = rs * Fmin_z[nl, :, k] * Po[nl]
        
        piLtot = piLo1 + piLo2 + piLo3 + piLo4
        LoF_[:, k] = piLtot / np.pi

    Fhem_ = np.sum(Fplu_z[0, :, :], axis=1)

    # --- Output (Additive) ---
    
    # Add Z modification to existing radiances
    rad['Lo_'][iwlfi_idx] += np.sum(LoF_, axis=1)
    
    # Add Z modification to reflectances/fluxes
    rad['rso'][iwlfi_idx] += LoF_[:, 0] / rad['Esun_'][iwlfi_idx]
    rad['rdo'][iwlfi_idx] += LoF_[:, 1] / rad['Esky_'][iwlfi_idx]
    
    rad['refl'][iwlfi_idx] = np.pi * rad['Lo_'][iwlfi_idx] / (rad['Esky_'][iwlfi_idx] + rad['Esun_'][iwlfi_idx])
    
    rad['Eout_'][iwlfi_idx] += Fhem_

    return rad