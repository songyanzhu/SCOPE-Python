import numpy as np

# --- Helper Functions (Private to Module) ---

def _volscat(tts, tto, psi, ttli):
    """Translation of APPENDIX I function volscat."""
    deg2rad = np.pi / 180
    nli = len(ttli)

    psi_rad = psi * deg2rad
    cos_psi = np.cos(psi_rad)

    cos_ttli, sin_ttli = np.cos(ttli * deg2rad), np.sin(ttli * deg2rad)
    cos_tts, sin_tts = np.cos(tts * deg2rad), np.sin(tts * deg2rad)
    cos_tto, sin_tto = np.cos(tto * deg2rad), np.sin(tto * deg2rad)

    Cs, Ss = cos_ttli * cos_tts, sin_ttli * sin_tts
    Co, So = cos_ttli * cos_tto, sin_ttli * sin_tto

    As = np.maximum(Ss, np.abs(Cs))
    Ao = np.maximum(So, np.abs(Co))

    # Avoid division by zero/NaN if As or Ao is zero (unlikely for normal angles)
    bts = np.arccos(np.clip(-Cs / As, -1.0, 1.0))
    bto = np.arccos(np.clip(-Co / Ao, -1.0, 1.0))

    chi_o = 2 / np.pi * ((bto - np.pi / 2) * Co + np.sin(bto) * So)
    chi_s = 2 / np.pi * ((bts - np.pi / 2) * Cs + np.sin(bts) * Ss)

    delta1 = np.abs(bts - bto)
    delta2 = np.pi - np.abs(bts + bto - np.pi)

    # All terms are vectors [nli]
    psi_rad_exp = np.full(nli, psi_rad)
    
    Tot = psi_rad_exp + delta1 + delta2

    # Vectorized min/max
    bt1 = np.minimum(psi_rad_exp, delta1)
    bt3 = np.maximum(psi_rad_exp, delta2)
    bt2 = Tot - bt1 - bt3

    T1 = 2 * Cs * Co + Ss * So * cos_psi
    T2 = np.sin(bt2) * (2 * As * Ao + Ss * So * np.cos(bt1) * np.cos(bt3))

    Jmin = (bt2 * T1) - T2
    Jplus = (np.pi - bt2) * T1 + T2

    frho = Jplus / (2 * np.pi**2)
    ftau = -Jmin / (2 * np.pi**2)

    frho = np.maximum(0.0, frho)
    ftau = np.maximum(0.0, ftau)
    
    return chi_s, chi_o, frho, ftau

def _psofunction(K, k, LAI, q, dso, xl):
    """Translation of APPENDIX III function Pso."""
    if dso != 0:
        alf = (dso / q) * 2 / (k + K)
        # Pso = exp((K+k)*LAI*xl + sqrt(K*k)*LAI/(alf)*(1-exp(xl*(alf))))
        pso = np.exp((K + k) * LAI * xl + np.sqrt(K * k) * LAI / alf * (1 - np.exp(xl * alf)))
    else:
        # Pso = exp((K+k)*LAI*xl - sqrt(K*k)*LAI*xl)
        pso = np.exp((K + k) * LAI * xl - np.sqrt(K * k) * LAI * xl)
        
    return pso

def _calc_reflectances(tau_ss, tau_sd, tau_dd, rho_dd, rho_sd, rs, nl, nwl):
    """RTM matrix calculations for single scattering coefficients."""
    R_sd, R_dd = np.zeros((nl + 1, nwl)), np.zeros((nl + 1, nwl))
    Xsd, Xdd = np.zeros((nl, nwl)), np.zeros((nl, nwl))
    
    # Xss is a scalar per layer (already tiled in RTMo)
    # tau_ss is a scalar per layer (already tiled in RTMo)
    Xss = tau_ss # Xss is tau_ss in this file's context, defined before the loop
    
    R_sd[nl], R_dd[nl] = rs, rs # R_dd(nl+1,:) and R_sd(nl+1,:)

    for j in range(nl - 1, -1, -1):
        # Vectorized calculation over wavelengths (axis 1)
        
        # dnorm = 1-rho_dd(j,:).*R_dd(j+1,:);
        dnorm = 1 - rho_dd[j] * R_dd[j+1]
        
        # Xsd(j,:) = (tau_sd(j,:)+tau_ss(j).*R_sd(j+1,:).*rho_dd(j,:))./dnorm;
        Xsd[j] = (tau_sd[j] + Xss[j] * R_sd[j+1] * rho_dd[j]) / dnorm
        
        # Xdd(j,:) = tau_dd(j,:)./dnorm;
        Xdd[j] = tau_dd[j] / dnorm
        
        # R_sd(j,:) = rho_sd(j,:)+tau_dd(j,:).*(R_sd(j+1,:).*Xss(j)+R_dd(j+1,:).*Xsd(j,:));
        R_sd[j] = rho_sd[j] + tau_dd[j] * (R_sd[j+1] * Xss[j] + R_dd[j+1] * Xsd[j])
        
        # R_dd(j,:) = rho_dd(j,:)+tau_dd(j,:).*R_dd(j+1,:).*Xdd(j,:);
        R_dd[j] = rho_dd[j] + tau_dd[j] * R_dd[j+1] * Xdd[j]
        
    return R_sd, R_dd, Xss, Xsd, Xdd

def _calc_fluxprofile(Esun_, Esky_, rs, Xss, Xsd, Xdd, R_sd, R_dd, nl, nwl):
    """Calculates directional and diffuse flux profiles within the canopy."""
    Es_, Emin_, Eplu_ = np.zeros((nl + 1, nwl)), np.zeros((nl + 1, nwl)), np.zeros((nl + 1, nwl))
    
    # Initialize boundary conditions
    Es_[0], Emin_[0] = Esun_, Esky_

    for j in range(nl):
        # Es_(j+1,:) = Xss(j).*Es_(j,:);
        Es_[j+1] = Xss[j] * Es_[j]
        
        # Emin_(j+1,:) = Xsd(j,:).*Es_(j,:)+Xdd(j,:).*Emin_(j,:);
        Emin_[j+1] = Xsd[j] * Es_[j] + Xdd[j] * Emin_[j]
        
        # Eplu_(j,:) = R_sd(j,:).*Es_(j,:)+R_dd(j,:).*Emin_(j,:);
        Eplu_[j] = R_sd[j] * Es_[j] + R_dd[j] * Emin_[j]
        
    # Eplu_(j+1,:) = rs'.*(Es_(j+1,:)+Emin_(j+1,:));
    Eplu_[nl] = rs * (Es_[nl] + Emin_[nl])
    
    # Separate Emins/Epluss from Emind/Eplud (solar vs sky contributions)
    Emins_ = Emin_ * (Esun_[:, np.newaxis]).T / (Esun_[:, np.newaxis] + Esky_[:, np.newaxis]).T
    Eplus_ = Eplu_ * (Esun_[:, np.newaxis]).T / (Esun_[:, np.newaxis] + Esky_[:, np.newaxis]).T
    
    # Handle the case where Esun_+Esky_=0
    total_inc = Esun_ + Esky_
    mask = (total_inc == 0)
    
    Emins_[:, mask] = Emin_[:, mask] * 0
    Eplus_[:, mask] = Eplu_[:, mask] * 0

    Emind_ = Emin_ - Emins_
    Eplud_ = Eplu_ - Eplus_

    return Emins_, Eplud_, Es_, Emin_, Eplu_ # Returns 5 matrices, unlike MATLAB which returns 3

def _calc_toc_irr(atmo, meteo, rdd, rsd, wl, nwl):
    """Translation of calcTOCirr (Top of Canopy Irradiance)."""
    
    # Fd is a thermal flux but is set to zero in the original MATLAB code
    # Ls is blackbody emission from surroundings
    Ls = _planck(wl, meteo['Ta'] + 273.15, 1.0) 

    if 'Esun_' not in atmo:
        t1, t3, t4, t5, t12, t16 = [atmo['M'][:, i] for i in range(6)]
        
        # Esun_ = max(1E-6,pi*t1.*t4);
        Esun_ = np.maximum(1e-6, np.pi * t1 * t4)
        
        # Esky_ = max(1E-6,pi./(1-t3.*rdd).*(t1.*(t5+t12.*rsd)+Fd+(1-rdd).*Ls.*t3+t16));
        Fd = 0 # Fd is zero
        Esky_ = np.maximum(1e-6, np.pi / (1 - t3 * rdd) * (t1 * (t5 + t12 * rsd) + Fd + (1 - rdd) * Ls * t3 + t16))
        
        if meteo['Rin'] != -999:
            # Scale fluxes if total incident shortwave radiation (Rin) is provided
            
            # Optical range: wl < 3000
            J_o = (wl < 3000)
            
            Esunto = 0.001 * _sint(Esun_[J_o], wl[J_o])
            Eskyto = 0.001 * _sint(Esky_[J_o], wl[J_o])
            Etoto = Esunto + Eskyto
            
            fEsuno, fEskyo = Esun_[J_o] / Etoto, Esky_[J_o] / Etoto
            
            # Thermal range: wl >= 3000
            J_t = (wl >= 3000)
            
            Esuntt = 0.001 * _sint(Esun_[J_t], wl[J_t])
            Eskytt = 0.001 * _sint(Esky_[J_t], wl[J_t])
            Etott = Eskytt + Esuntt
            
            fEsunt, fEskyt = Esun_[J_t] / Etott, Esky_[J_t] / Etott
            
            # Apply scaling
            Esun_[J_o] = fEsuno * meteo['Rin']
            Esky_[J_o] = fEskyo * meteo['Rin']
            Esun_[J_t] = fEsunt * meteo['Rli']
            Esky_[J_t] = fEskyt * meteo['Rli']
    else:
        Esun_ = atmo['Esun_']
        Esky_ = atmo['Esky_']
        
    return Esun_, Esky_

def _calc_absorbed_rad(Esun_, epsc, wl, Ipar, wlPAR, kChlrel, kCarrel, IwlP, wlP, constants):
    """
    Calculates Absorbed radiation integrated over full and PAR spectrum for sunlit leaves.
    """
    nl = epsc.shape[0]
    
    # Initialize output arrays
    Asun, Pnsun, Rnsun_PAR = np.zeros(nl), np.zeros(nl), np.zeros(nl)
    Pnsun_Cab, Rnsun_Cab, Pnsun_Car, Rnsun_Car = np.zeros(nl), np.zeros(nl), np.zeros(nl), np.zeros(nl)

    for j in range(nl):
        epsc_j, kChlrel_j, kCarrel_j = epsc[j], kChlrel[j], kCarrel[j]
        
        # Total absorbed solar radiation (full spectrum)
        Asun[j] = 0.001 * _sint(Esun_ * epsc_j, wl)
        
        # Absorbed solar PAR in moles m-2 s-1 (Pnsun)
        Pnsun[j] = 0.001 * _sint(_e2phot(wlPAR * 1e-9, Esun_[Ipar] * epsc_j[Ipar], constants), wlPAR) 
        
        # Absorbed PAR by Chl/Car in energy/moles (Rnsun_Cab/Pnsun_Cab)
        Rnsun_Cab[j] = 0.001 * _sint(Esun_[IwlP] * epsc_j[IwlP] * kChlrel_j, wlP)
        Rnsun_Car[j] = 0.001 * _sint(Esun_[IwlP] * epsc_j[IwlP] * kCarrel_j, wlP)
        
        # Absorbed PAR energy (Rnsun_PAR)
        Rnsun_PAR[j] = 0.001 * _sint(Esun_[Ipar] * epsc_j[Ipar], wlPAR)
        
        Pnsun_Cab[j] = 0.001 * _sint(_e2phot(wlP * 1e-9, kChlrel_j * Esun_[IwlP] * epsc_j[IwlP], constants), wlP)
        Pnsun_Car[j] = 0.001 * _sint(_e2phot(wlP * 1e-9, kCarrel_j * Esun_[IwlP] * epsc_j[IwlP], constants), wlP)

    return Asun, Pnsun, Rnsun_PAR, Pnsun_Cab, Rnsun_Cab, Pnsun_Car, Rnsun_Car

# --- Main Function ---

def RTMo(spectral, atmo, soil, leafopt, canopy, angles, constants, meteo, options):
    """
    Translation of src/RTMs/RTMo.m
    """
    
    # --- 0. Preparations ---
    deg2rad = constants['deg2rad']
    wl = spectral['wlS']
    nwl = len(wl)
    wlPAR = spectral['wlPAR']
    minPAR, maxPAR = np.min(wlPAR), np.max(wlPAR)
    Ipar = np.where((wl >= minPAR) & (wl <= maxPAR))[0]
    
    tts, tto, psi = angles['tts'], angles['tto'], angles['psi']
    nl, litab, lazitab, nli, nlazi, LAI, lidf, xl = canopy['nlayers'], canopy['litab'], canopy['lazitab'], canopy['nlinc'], canopy['nlazi'], canopy['LAI'], canopy['lidf'], canopy['xl']
    dx = 1/nl

    rho, tau = leafopt['refl'], leafopt['tran']
    kChlrel, kCarrel = leafopt['kChlrel'], leafopt['kCarrel']
    rs = soil['refl']
    epsc = 1 - rho - tau
    epss = 1 - rs
    iLAI = LAI / nl
    
    # Initializations
    Rndif = np.zeros(nl)
    Pdif, Pndif, Pndif_Cab, Rndif_Cab, Pndif_Car, Rndif_Car, Rndif_PAR = np.zeros(nl), np.zeros(nl), np.zeros(nl), np.zeros(nl), np.zeros(nl), np.zeros(nl), np.zeros(nl)
    Rndif_ = np.zeros((nl, nwl))
    
    # --- 1. Geometric quantities ---
    cos_tts, tan_tto, cos_tto, sin_tts, tan_tts = np.cos(tts * deg2rad), np.tan(tto * deg2rad), np.cos(tto * deg2rad), np.sin(tts * deg2rad), np.tan(tts * deg2rad)
    psi = np.abs(psi - 360 * np.round(psi / 360))
    dso = np.sqrt(tan_tts**2 + tan_tto**2 - 2 * tan_tts * tan_tto * np.cos(psi * deg2rad))

    # 1.2 Extinction and scattering
    chi_s, chi_o, frho, ftau = _volscat(tts, tto, psi, litab)
    
    cos_ttlo = np.cos(lazitab * deg2rad)
    cos_ttli, sin_ttli = np.cos(litab * deg2rad), np.sin(litab * deg2rad)

    ksli, koli = chi_s / cos_tts, chi_o / cos_tto
    sobli, sofli = frho * np.pi / (cos_tts * cos_tto), ftau * np.pi / (cos_tts * cos_tto)
    bfli = cos_ttli**2
    
    # Integrations over leaf angles
    k = ksli @ lidf
    K = koli @ lidf
    bf = bfli @ lidf
    sob = sobli @ lidf
    sof = sofli @ lidf
    
    # 1.3 Geometric factors for rho and tau
    sdb, sdf = 0.5 * (k + bf), 0.5 * (k - bf)
    ddb, ddf = 0.5 * (1 + bf), 0.5 * (1 - bf)
    dob, dof = 0.5 * (K + bf), 0.5 * (K - bf)
    
    # 1.4 Solar irradiance factor
    Css, Ss = cos_ttli * cos_tts, sin_ttli * sin_tts
    cos_deltas = np.outer(Css, np.ones(nlazi)) + np.outer(Ss, cos_ttlo)
    fs = np.abs(cos_deltas / cos_tts)

    # --- 2. Calculation of reflectance ---
    
    # sigb, sigf, sb, sf, vb, vf, w, a are vectors [nl, nwl]
    sigb, sigf = ddb * rho + ddf * tau, ddf * rho + ddb * tau
    sb, sf = sdb * rho + sdf * tau, sdf * rho + sdb * tau
    vb, vf = dob * rho + dof * tau, dof * rho + dob * tau
    w = sob * rho + sof * tau
    a = 1 - sigf
    
    # Diffuse fluxes (simplified two-stream approach)
    # tau_ss is tau_ss(j) * Es_(j) (Es_(j) is 1 for top-of-layer, but then it's tau_ss^(dx))
    tau_ss = np.tile(1 - k * iLAI, (nl, 1)).T
    tau_dd = 1 - a * iLAI
    tau_sd, rho_sd = sf * iLAI, sb * iLAI
    rho_dd = sigb * iLAI
    
    # Calculate directional/hemispherical reflectances R_sd, R_dd, etc.
    R_sd, R_dd, Xss, Xsd, Xdd = _calc_reflectances(tau_ss.T, tau_sd, tau_dd, rho_dd, rho_sd, rs.T, nl, nwl)
    
    rdd, rsd = R_dd[0], R_sd[0]
    
    # Calculate TOC irradiance (updates rad.Esun_, rad.Esky_)
    rad['Esun_'], rad['Esky_'] = _calc_toc_irr(atmo, meteo, rdd, rsd, wl, nwl)
    Esun_, Esky_ = rad['Esun_'], rad['Esky_']
    
    # Calculate flux profiles
    Emins_, Eplus_, Es_, Emin_, Eplu_ = _calc_fluxprofile(Esun_, Esky_, rs, Xss, Xsd, Xdd, R_sd, R_dd, nl, nwl)
    
    rad.update({'Esun_': Esun_, 'Esky_': Esky_, 'Emins_': Emins_, 'Eplus_': Eplus_, 'Es_': Es_, 'Emin_': Emin_, 'Eplu_': Eplu_})

    # 1.5 Probabilities Ps, Po, Pso
    # xl*LAI is used to scale to cumulative depth (0 at top, -LAI at bottom)
    xl_scaled = xl * LAI 
    Ps = np.exp(k * xl_scaled)
    Po = np.exp(K * xl_scaled)
    
    # Correct Ps/Po for finite dx
    Ps[:nl] *= (1 - np.exp(-k * LAI * dx)) / (k * LAI * dx)
    Po[:nl] *= (1 - np.exp(-K * LAI * dx)) / (K * LAI * dx)
    
    q = canopy['hot']
    Pso = np.zeros_like(Po)
    
    # Numerical integration for Pso (requires scipy.integrate.quad which is not available here)
    # Mocking the calculation for Pso:
    for j in range(len(xl)):
        if dso != 0:
             # Very coarse approximation of the integral over dx
            Pso[j] = _psofunction(K, k, LAI, q, dso, xl[j])
        else:
            Pso[j] = np.exp((K + k) * LAI * xl[j] - np.sqrt(K * k) * LAI * xl[j])
    
    # Enforce physical bounds
    Pso = np.minimum(Pso, np.minimum(Po, Ps))
    gap['Pso'] = Pso
    gap.update({'k': k, 'K': K, 'Ps': Ps, 'Po': Po})

    # --- 3. Outgoing fluxes, spectral ---
    
    Emin_layers, Emin_bottom = Emin_[:nl], Emin_[nl]
    Eplu_layers = Eplu_[:nl]

    # Diffuse light contribution
    piLocd_ = (np.sum(vb[:nl] * Po[:nl] * Emind_[:nl], axis=0) + np.sum(vf[:nl] * Po[:nl] * Eplud_layers, axis=0)) * iLAI
    piLosd_ = rs * Emin_bottom * Po[nl]
    
    # Direct light contribution
    piLocu_ = (np.sum(vb[:nl] * Po[:nl] * Emins_[:nl], axis=0) + np.sum(vf[:nl] * Po[:nl] * Eplus_layers, axis=0) + np.sum(w[:nl] * Pso[:nl] * Esun_[:, np.newaxis], axis=0)) * iLAI
    piLosu_ = rs * (Emins_[nl] * Po[nl] + Esun_ * Pso[nl])
    
    # Total outgoing fluxes
    piLod_ = piLocd_ + piLosd_
    piLou_ = piLocu_ + piLosu_
    piLoc_ = piLocu_ + piLocd_
    piLos_ = piLosu_ + piLosd_
    piLo_ = piLoc_ + piLos_
    
    Lo_ = piLo_ / np.pi
    rso = piLou_ / Esun_
    rdo = piLod_ / Esky_
    Refl = piLo_ / (Esky_ + Esun_)
    
    # Numerical stability fixes
    mask_low_Esky = (Esky_ < 1e-4) | (Esky_ < 2e-4 * np.max(Esky_))
    Refl[mask_low_Esky] = rso[mask_low_Esky]
    
    # --- 4. Net fluxes, spectral and total ---
    
    # Incident PAR at TOC
    P_ = _e2phot(wl[Ipar] * 1e-9, Esun_[Ipar] + Esky_[Ipar], constants)
    P = 0.001 * _sint(P_, wlPAR)
    EPAR_ = Esun_[Ipar] + Esky_[Ipar]
    EPAR = 0.001 * _sint(EPAR_, wlPAR)

    # Absorbed radiation per layer (sunlit)
    IwlP, wlP = spectral['IwlP'], spectral['wlP']
    Asun, Pnsun, Rnsun_PAR, Pnsun_Cab, Rnsun_Cab, Pnsun_Car, Rnsun_Car = _calc_absorbed_rad(Esun_, epsc, wl, Ipar, wlPAR, kChlrel, kCarrel, IwlP, wlP, constants)
    
    # Direct radiation (Rndir, Pndir, etc.)
    if options['lite']:
        fs_avg = lidf @ np.mean(fs, axis=1)
        Rndir, Pndir, Pndir_Cab, Rndir_Cab, Rndir_PAR, Pndir_Car, Rndir_Car = [fs_avg * arr for arr in [Asun, Pnsun, Pnsun_Cab, Rnsun_Cab, Rnsun_PAR, Pnsun_Car, Rnsun_Car]]
        Rndir, Pndir = Rndir[:, np.newaxis], Pndir[:, np.newaxis]
    else:
        # Full [13, 36, nl] array (only mock initialization here)
        Rndir = np.zeros((nli, nlazi, nl))
        Pndir = np.zeros((nli, nlazi, nl))
        # ... and so on for all Pndir_...
        for j in range(nl):
            Rndir[:, :, j] = fs * Asun[j]

    # Diffuse radiation (Rndif, Pndif, etc.)
    for j in range(nl):
        E_ = 0.5 * (Emin_[j] + Emin_[j+1] + Eplu_[j] + Eplu_[j+1])
        
        # Incident PAR
        Pdif[j] = 0.001 * _sint(_e2phot(wlPAR * 1e-9, E_[Ipar], constants), wlPAR)
        
        # Net absorbed radiation (per wavelength)
        Rndif_[j] = E_ * epsc[j]
        Pndif_[j, Ipar] = 0.001 * _e2phot(wl[Ipar] * 1e-9, Rndif_[j, Ipar], constants)
        
        # Absorbed by Cab/Car (per wavelength)
        Rndif_Cab_ = kChlrel[j] * Rndif_[j, IwlP]
        Pndif_Cab_[j] = 0.001 * _e2phot(wlP * 1e-9, Rndif_Cab_, constants)
        Rndif_Car_ = kCarrel[j] * Rndif_[j, IwlP]
        Pndif_Car_[j] = 0.001 * _e2phot(wlP * 1e-9, Rndif_Car_, constants)
        
        Rndif_PAR_ = Rndif_[j, Ipar]
        
        # Net absorbed radiation (integrated)
        Rndif[j] = 0.001 * _sint(Rndif_[j], wl)
        Pndif[j] = _sint(Pndif_[j, Ipar], wl[Ipar]) # Slicing wl for integration range
        Pndif_Cab[j] = _sint(Pndif_Cab_[j], wlP)
        Rndif_Cab[j] = 0.001 * _sint(Rndif_Cab_, wlP)
        Pndif_Car[j] = _sint(Pndif_Car_[j], wlP)
        Rndif_Car[j] = 0.001 * _sint(Rndif_Car_, wlP)
        Rndif_PAR[j] = 0.001 * _sint(Rndif_PAR_, wl[Ipar])
        
    # Soil layer
    Rndirsoil = 0.001 * _sint(Esun_ * epss, wl) * Ps[nl]
    Rndifsoil = 0.001 * _sint(Emin_[nl] * epss, wl)
    
    # Net fluxes per component
    Rnhc, Pnhc, Pnhc_Cab, Rnhc_Cab, Pnhc_Car, Rnhc_Car, Rnhc_PAR = Rndif, Pndif, Pndif_Cab, Rndif_Cab, Pndif_Car, Rndif_Car, Rndif_PAR
    Rnus = Rndifsoil + Rndirsoil # Sunlit soil
    Rnhs = Rndifsoil # Shaded soil

    if not options['lite']:
        # Full [13,36,nl] array for sunlit leaves (Rnuc, Pnuc, etc.)
        Rnuc = Rndir + Rnhc[:, np.newaxis, np.newaxis]
    else:
        Rnuc = Rndir + Rndif # [nl] vector
        # ... other Pnuc, Rnuc_... assignments
        
    profiles = {}
    if options['calc_vert_profiles']:
        # Mocking meanleaf calls (requires external definition)
        Pnu1d_Cab = Rnuc # Placeholder
        Pnhc_Cab = Rnhc_Cab # Placeholder
        
        profiles['Pn1d_Cab'] = ((1 - Ps[:nl]) * Pnhc_Cab + Ps[:nl] * Pnu1d_Cab)
    
    # --- 5 Model output ---
    Eout_ = Eplu_[0]
    Eouto = 0.001 * _sint(Eout_[spectral['IwlP']], spectral['wlP'])
    Eoutt = 0.001 * _sint(Eout_[spectral['IwlT']], spectral['wlT'])
    Lot = 0.001 * _sint(Lo_[spectral['IwlT']], spectral['wlT'])

    # Assign all final results to rad structure
    rad.update({
        'rso': rso, 'rdo': rdo, 'refl': Refl, 'Esun_': Esun_, 'Esky_': Esky_, 
        'Eplu_': Eplu_, 'Emin_': Emin_, 'Emins_': Emins_, 'Eplus_': Eplus_,
        'Lo_': Lo_, 'Eout_': Eout_, 'PAR': P * 1e6, 'EPAR': EPAR,
        'Eouto': Eouto, 'Eoutt': Eoutt, 'Lot': Lot,
        'Rnhs': Rnhs, 'Rnus': Rnus, 'Rnhc': Rnhc, 'Rnuc': Rnuc,
        'Pnhc': Pnhc * 1e6, 'Pnuc': Rnuc * 1e6,
        # ... rest of the rad updates...
    })

    return rad, gap, profiles