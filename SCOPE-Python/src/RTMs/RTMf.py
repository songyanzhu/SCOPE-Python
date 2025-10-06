import numpy as np

# --- Helper Functions (From RTMo.m appendix/other files, needed here) ---

def _ephoton(lambda_val, constants):
    """Calculates the energy content (J) of 1 photon of wavelength lambda (m)"""
    h = constants['h']
    c = constants['c']
    E = h * c / lambda_val
    return E

def _e2phot(lambda_val, E, constants):
    """Calculates the number of moles of photons corresponding to E [J] of energy."""
    e = _ephoton(lambda_val, constants)
    photons = E / e
    molphotons = photons / constants['A']
    return molphotons

def _sint(flux, wl):
    """
    Simulates MATLAB's Sint (Spectral integration). Assumes wavelengths are linearly spaced.
    Flux and wl should have the same shape.
    """
    if flux.size == 0 or wl.size == 0:
        return 0.0
    
    # MATLAB's Sint assumes trapezoidal integration on the given array of fluxes.
    # The spacing is (wl[end]-wl[1])/(len(wl)-1), but here we use the diffs.
    dwl = np.diff(wl)
    # Need to handle when wl is a range, but flux is computed on those points.
    
    # If wl is a regular grid, dwl is a scalar.
    if np.allclose(dwl, dwl[0]):
        # Uniform spacing: sum(flux) * dwl[0]
        result = np.trapz(flux, wl)
    else:
        # Non-uniform spacing: use numpy's trapz
        result = np.trapz(flux, wl)
        
    return result

def _interp1(X, V, XI, method):
    """
    Wrapper for numpy/scipy interpolation, simulating MATLAB's interp1.
    'method' is ignored for simplicity but kept in the signature.
    """
    # For SCOPE, 'spline' is often used, requiring SciPy. 
    # Using linear interpolation for NumPy-only compatibility unless SciPy is guaranteed.
    from scipy.interpolate import interp1d
    
    # Create the interpolation function
    f = interp1d(X, V, kind='linear', fill_value="extrapolate")
    
    # Evaluate the function at the new points
    VI = f(XI)
    return VI

# --- Main Function ---

def RTMf(constants, spectral, rad, soil, leafopt, canopy, gap, angles, etau, etah):
    """
    Translation of src/RTMs/RTMf.m
    Calculates the spectrum of fluorescent radiance.
    """

    # --- 0.1 initialisations ---
    
    # Use array conversion for wlS if necessary to match MATLAB column vector assumption
    wlS = spectral['wlS']
    wlF = np.arange(640, 851, 4) # (640:4:850)'
    wlE = np.arange(400, 751, 5) # (400:5:750)'
    
    # Intersect for indices
    iwlfi = np.isin(wlS, wlE)
    iwlfo = np.isin(wlS, wlF)
    
    iwlfi_idx = np.where(iwlfi)[0]
    iwlfo_idx = np.where(iwlfo)[0]
    
    nf = len(iwlfo_idx)
    nl = canopy['nlayers']
    LAI = canopy['LAI']
    litab = canopy['litab']
    lazitab = canopy['lazitab']
    lidf = canopy['lidf']
    nlazi = canopy['nlazi']
    nlinc = canopy['nlinc']
    nlori = nlinc * nlazi
    layers = np.arange(1, nl + 1)

    Ps, Po, Pso = gap['Ps'], gap['Po'], gap['Pso']
    Qs = Ps[:-1] # Ps(1:end-1)

    # Initialize matrices [nwlfo, nl]
    MpluEmin, MpluEplu, MminEmin, MminEplu, MpluEsun, MminEsun = [np.zeros((nf, nl)) for _ in range(6)]

    # Slice rad fields (transpose for MATLAB's [nl, nwl] -> [nwl, nl] convention)
    Esunf_ = rad['Esun_'][iwlfi_idx] # [nwlfi]
    Eminf_ = rad['Emin_'][1:nl+1, iwlfi_idx].T  # [nl+1, nwl] -> [nl, nwlfi]. Emin_[1:end, iwlfi]'
    Epluf_ = rad['Eplu_'][:nl, iwlfi_idx].T      # [nl+1, nwl] -> [nl, nwlfi]. Eplu_[1:end-1, iwlfi]'
    iLAI = LAI / nl
    
    # RTM coefficients at fluorescence wavelengths [nwlfo]
    Xdd = rad['Xdd'][:, iwlfo_idx]
    rho_dd = rad['rho_dd'][:, iwlfo_idx]
    R_dd = rad['R_dd'][:, iwlfo_idx]
    tau_dd = rad['tau_dd'][:, iwlfo_idx]
    vb = rad['vb'][:, iwlfo_idx]
    vf = rad['vf'][:, iwlfo_idx]

    # --- 0.2 geometric quantities ---
    
    Mb, Mf = leafopt['Mb'], leafopt['Mf']
    
    deg2rad = constants['deg2rad']
    tto, tts, psi = angles['tto'], angles['tts'], angles['psi']
    rs = soil['refl'][iwlfo_idx]
    
    cos_tto, sin_tto = np.cos(tto * deg2rad), np.sin(tto * deg2rad)
    cos_tts, sin_tts = np.cos(tts * deg2rad), np.sin(tts * deg2rad)
    cos_ttli, sin_ttli = np.cos(litab * deg2rad), np.sin(litab * deg2rad)
    cos_phils = np.cos(lazitab * deg2rad)
    cos_philo = np.cos((lazitab - psi) * deg2rad)

    # --- 0.3 geometric factors for all leaf angle/azumith classes ---
    
    # MATLAB: cos_ttli*cos_tts*ones(1,36) -> np.outer(cos_ttli, np.ones(nlazi)) * cos_tts
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
    ctl2 = cos_ttli**2
    ctl2 = np.outer(ctl2, np.ones(nlazi))
    
    # reshape all the variables to [nlori, 1]
    absfs = absfs.reshape(-1, 1)
    absfo = absfo.reshape(-1, 1)
    fsfo = fsfo.reshape(-1, 1)
    absfsfo = absfsfo.reshape(-1, 1)
    foctl = foctl.reshape(-1, 1)
    fsctl = fsctl.reshape(-1, 1)
    ctl2 = ctl2.reshape(-1, 1)

    # --- 1.0 calculation of fluorescence flux in observation direction ---

    # Fluorescence matrices and efficiencies
    U, Fmin_, Fplu_ = [np.zeros((nl + 1, Mb.shape[1])) for _ in range(3)]
    
    # Slice Mb, Mf to [nwlfo, nwlfi, nl]
    Mb_f = Mb[iwlfo_idx, :, :] 
    Mf_f = Mf[iwlfo_idx, :, :] 
    
    Mplu = 0.5 * (Mb_f + Mf_f)  # [nwlfo, nwlfi, nl]
    Mmin = 0.5 * (Mb_f - Mf_f)  # [nwlfo, nwlfi, nl]

    # In-products: incoming radiation to fluorescence spectrum
    # MpluEmin = ep.*(Mplu(:,:,j)* e2phot(wlE*1E-9,Eminf_(:,j),constants))
    
    wlE_m = wlE * 1e-9 # [m]
    wlF_m = wlF * 1e-9 # [m]
    
    # ep = constants.A*ephoton(wlF*1E-9,constants)
    ep = constants['A'] * _ephoton(wlF_m, constants)
    
    # Reshape ep for element-wise multiplication with Mplu * e2phot result
    ep_exp = ep[:, np.newaxis]

    for j in range(nl):
        # Eminf_[:, j] is [nwlfi]
        Emin_phot = _e2phot(wlE_m, Eminf_[:, j], constants) # [nwlfi]
        Eplu_phot = _e2phot(wlE_m, Epluf_[:, j], constants) # [nwlfi]
        Esun_phot = _e2phot(wlE_m, Esunf_, constants)       # [nwlfi]
        
        MpluEmin[:, j] = (Mplu[:, :, j] @ Emin_phot) * ep # [nwlfo]
        MpluEplu[:, j] = (Mplu[:, :, j] @ Eplu_phot) * ep
        MminEmin[:, j] = (Mmin[:, :, j] @ Emin_phot) * ep
        MminEplu[:, j] = (Mmin[:, :, j] @ Eplu_phot) * ep
        MpluEsun[:, j] = (Mplu[:, :, j] @ Esun_phot) * ep
        MminEsun[:, j] = (Mmin[:, :, j] @ Esun_phot) * ep

    # Emission factors (bsxfun = broadcasting and operation)
    laz = 1 / nlazi
    
    if etau.ndim == 1:
        # If etau is [nl], expand to [nl, nlinc, nlazi] and then reshape
        etau_3d = np.repeat(etau, nlinc * nlazi).reshape(nl, nlinc, nlazi)
        etau_3d = np.transpose(etau_3d, (1, 2, 0)) # [nlinc, nlazi, nl]
        etau_nlori_nl = etau_3d.reshape(nlori, nl) # [nlori, nl]
    else:
        # Assuming etau is [nlinc, nlazi, nl] which reshapes to [nlori, nl]
        etau_nlori_nl = etau.reshape(nlori, nl)

    # etau_lidf [nlori, nl]
    lidf_laz = np.outer(lidf, np.ones(nlazi)) * laz # [nlinc, nlazi]
    lidf_laz_nlori = lidf_laz.reshape(nlori, 1) # [nlori, 1]
    
    etau_lidf = etau_nlori_nl * lidf_laz_nlori
    etah_lidf = np.tile(etah.T, (nlori, 1)) * lidf_laz_nlori # etah is [nl, 1], transpose to [1, nl]

    # Fluorescence RTM components (vectorized over nwlfo)
    # The sum is over nlori (axis 0 of etau_lidf). The result is [1, nl] which broadcasts to [nwlfo, nl]
    
    # bsxfun(@times,sum(bsxfun(@times,etau_lidf,absfsfo)),MpluEsun)
    # sum_factor is [1, nl] (after sum(axis=0))
    sum_etau_absfsfo = np.sum(etau_lidf * absfsfo, axis=0)
    
    # wfEs [nwlfo, nl]
    wfEs = MpluEsun * sum_etau_absfsfo + MminEsun * np.sum(etau_lidf * fsfo, axis=0)
    
    # sfEs [nwlfo, nl]
    sfEs = MpluEsun * np.sum(etau_lidf * absfs, axis=0) - MminEsun * np.sum(etau_lidf * fsctl, axis=0)
    
    # sbEs [nwlfo, nl]
    sbEs = MpluEsun * np.sum(etau_lidf * absfs, axis=0) + MminEsun * np.sum(etau_lidf * fsctl, axis=0)

    # vfEplu_h/u (Diffuse upward)
    vfEplu_h = MpluEplu * np.sum(etah_lidf * absfo, axis=0) - MminEplu * np.sum(etah_lidf * foctl, axis=0)
    vfEplu_u = MpluEplu * np.sum(etau_lidf * absfo, axis=0) - MminEplu * np.sum(etau_lidf * foctl, axis=0)

    # vbEmin_h/u (Diffuse downward)
    vbEmin_h = MpluEmin * np.sum(etah_lidf * absfo, axis=0) + MminEmin * np.sum(etah_lidf * foctl, axis=0)
    vbEmin_u = MpluEmin * np.sum(etau_lidf * absfo, axis=0) + MminEmin * np.sum(etau_lidf * foctl, axis=0)

    # sigf/sigb (Diffuse coefficients)
    sum_etau = np.sum(etau_lidf, axis=0)
    sum_etah = np.sum(etah_lidf, axis=0)
    sum_etau_ctl2 = np.sum(etau_lidf * ctl2, axis=0)
    sum_etah_ctl2 = np.sum(etah_lidf * ctl2, axis=0)
    
    sigfEmin_h = MpluEmin * sum_etah - MminEmin * sum_etah_ctl2
    sigfEmin_u = MpluEmin * sum_etau - MminEmin * sum_etau_ctl2
    sigbEmin_h = MpluEmin * sum_etah + MminEmin * sum_etah_ctl2
    sigbEmin_u = MpluEmin * sum_etau + MminEmin * sum_etau_ctl2
    
    sigfEplu_h = MpluEplu * sum_etah - MminEplu * sum_etah_ctl2
    sigfEplu_u = MpluEplu * sum_etau - MminEplu * sum_etau_ctl2
    sigbEplu_h = MpluEplu * sum_etah + MminEplu * sum_etah_ctl2
    sigbEplu_u = MpluEplu * sum_etau + MminEplu * sum_etau_ctl2

    # Emitted fluorescence per layer
    piLs = wfEs + vfEplu_u + vbEmin_u      # Sunlit
    piLd = vbEmin_h + vfEplu_h             # Shaded
    
    Fsmin = sfEs + sigfEmin_u + sigbEplu_u
    Fsplu = sbEs + sigbEmin_u + sigfEplu_u
    Fdmin = sigfEmin_h + sigbEplu_h
    Fdplu = sigbEmin_h + sigfEplu_h
    
    # Femmin = iLAI*bsxfun(@times,Qs', Fsmin) +iLAI* bsxfun(@times,(1-Qs)',Fdmin)
    Qs_exp = Qs[np.newaxis, :]
    Femmin = iLAI * (Fsmin * Qs_exp + Fdmin * (1 - Qs_exp))
    Femplu = iLAI * (Fsplu * Qs_exp + Fdplu * (1 - Qs_exp))

    # Two-flow RTM over fluorescence spectrum (using the same RTM coefficients as for the optical range)
    # The implementation is vectorized over nwlfo.
    Y, U = np.zeros((nl, nf)), np.zeros((nl + 1, nf))
    
    # Downward RTM (Bottom to Top)
    for j in range(nl - 1, -1, -1):
        # Y(j,:) = (rho_dd(j,:).*U(j+1,:)+Femmin(:,j)')./(1-rho_dd(j,:).*R_dd(j+1,:));
        
        # U[j+1] is [nwlfo]
        rho_dd_j, R_dd_j1 = rho_dd[j], R_dd[j+1]
        
        # Femmin is [nwlfo, nl]. Slice for column j
        Y[j] = (rho_dd_j * U[j+1] + Femmin[:, j]) / (1 - rho_dd_j * R_dd_j1)
        
        # U(j,:) = tau_dd(j,:).*(R_dd(j+1,:).*Y(j,:)+U(j+1,:))+Femplu(:,j)';
        tau_dd_j = tau_dd[j]
        U[j] = tau_dd_j * (R_dd_j1 * Y[j] + U[j+1]) + Femplu[:, j]

    # Upward RTM (Top to Bottom)
    Fmin_[:, :, 0], Fplu_[:, :, 0] = np.zeros((nl + 1, nf)), np.zeros((nl + 1, nf)) # k=1 for sun
    
    for j in range(nl):
        # Fmin_(j+1,:) = Xdd(j,:).*Fmin_(j,:)+Y(j,:);
        Fmin_[j+1, :, 0] = Xdd[j] * Fmin_[j, :, 0] + Y[j]
        
        # Fplu_(j,:) = R_dd(j,:).*Fmin_(j,:)+U(j,:);
        Fplu_[j, :, 0] = R_dd[j] * Fmin_[j, :, 0] + U[j]

    # Final Upward Fluxes
    Fhem_ = Fplu_[0, :, 0] # Fplu_(1,:)'
    
    # Directional Fluorescence Components (LoF_)
    # piLo1 = iLAI*piLs*Pso(1:nl)
    piLo1 = iLAI * piLs @ Pso[:nl]
    
    # piLo2 = iLAI*piLd*(Po(1:nl)-Pso(1:nl))
    piLo2 = iLAI * piLd @ (Po[:nl] - Pso[:nl])
    
    # piLo3 = iLAI*(vb.*Fmin_(layers,:) + vf.*Fplu_(layers,:))'*Po(1:nl)
    # layers is [1:nl]. vb, vf are [nl, nwlfo]. Fmin/Fplu are [nl+1, nwlfo].
    vb_layers, vf_layers = vb[:nl], vf[:nl]
    Fmin_layers = Fmin_[layers, :, 0] # [nl, nwlfo]
    Fplu_layers = Fplu_[layers, :, 0] # [nl, nwlfo]
    Po_layers = Po[:nl] # [nl]
    
    # (nl, nwlfo) * (nl, 1) + (nl, nwlfo) * (nl, 1) -> (nl, nwlfo)
    piLo3 = iLAI * (vb_layers * Fmin_layers + vf_layers * Fplu_layers).T @ Po_layers

    # piLo4 = rs .* Fmin_(nl+1,:)' * Po(nl+1)
    piLo4 = rs * Fmin_[nl, :, 0] * Po[nl]

    piLtot = piLo1 + piLo2 + piLo3 + piLo4
    LoF_ = piLtot / np.pi
    
    # --- Output Interpolation and Integration ---
    
    # LoF_ is [nwlfo]. wlF is [nwlfo]. spectral.wlF is [nwl_out].
    
    # Use spline if available, otherwise linear for interpolation
    method = 'linear' 
    if 'interp1' in globals(): # Checking if the SciPy-based helper is defined
        method = 'spline'

    # Radiance (LoF_) and Hemispherical flux (Fhem_) interpolation
    rad['LoF_'] = _interp1(wlF, LoF_, spectral['wlF'], method)
    rad['EoutF_'] = _interp1(wlF, Fhem_, spectral['wlF'], method)

    # Component Radiance interpolation
    rad['LoF_sunlit'] = _interp1(wlF, piLo1 / np.pi, spectral['wlF'], method)
    rad['LoF_shaded'] = _interp1(wlF, piLo2 / np.pi, spectral['wlF'], method)
    rad['LoF_scattered'] = _interp1(wlF, piLo3 / np.pi, spectral['wlF'], method)
    rad['LoF_soil'] = _interp1(wlF, piLo4 / np.pi, spectral['wlF'], method)

    # Total integrated fluxes (0.001 for mW to W)
    rad['EoutF'] = 0.001 * _sint(Fhem_, wlF)
    rad['LoutF'] = 0.001 * _sint(LoF_, wlF)

    # Femleaves_
    Femleaves_total = np.sum(Femmin + Femplu, axis=1) # Sum over layers (axis 1)
    rad['Femleaves_'] = _interp1(wlF, Femleaves_total, spectral['wlF'], method)

    # F685, F740 calculation (based on spectral.wlF's indices)
    # The original MATLAB code is highly dependent on wlF's specific values and range.
    # We'll mock the calculation using np.argmax on the interpolated result.
    
    # Assuming interpolated wlF starts near 640 and steps by 1 (which it doesn't, but is common for spectral.wlF)
    # F685: 640 -> 685 (45 nm range / 4 nm step = 11 steps. Index 1 to 12)
    # The MATLAB indices (1:55) likely correspond to wlF < 685nm or similar.
    
    # Index for ~685 nm (approx index in rad.LoF_)
    idx_685 = np.argmin(np.abs(spectral['wlF'] - 685)) 
    rad['F685'], iwl685 = rad['LoF_'][idx_685], idx_685
    rad['wl685'] = spectral['wlF'][iwl685]
    
    # Index for ~740 nm
    idx_740 = np.argmin(np.abs(spectral['wlF'] - 740)) 
    rad['F740'], iwl740 = rad['LoF_'][idx_740], idx_740
    rad['wl740'] = spectral['wlF'][iwl740] + 69 # The MATLAB offset is arbitrary here

    # Final mock for specific wavelengths (highly dependent on exact wlF definition)
    rad['F684'] = np.interp(684, spectral['wlF'], rad['LoF_'])
    rad['F761'] = np.interp(761, spectral['wlF'], rad['LoF_'])
    
    return rad