import numpy as np

def calczenithangle(Doy, t, Omega_g=210, Fi_gm=30, Long=13.75, Lat=45.5):
    """
    Calculates various solar and slope zenith and azimuth angles.
    
    :param Doy: Day of the year.
    :param t: Time of the day (hours, GMT).
    :param Omega_g: Slope azimuth angle (deg, default 210).
    :param Fi_gm: Slope of the surface (deg, default 30).
    :param Long: Longitude (decimal, default 13.75).
    :param Lat: Latitude (decimal, default 45.5).
    :return: Fi_s, Fi_gs, Fi_g, Omega_s (all in radians).
    """

    # convert angles into radials
    # G = (Doy-1)/365*2*pi
    G = (Doy - 1) / 365 * 2 * np.pi
    
    # Convert degrees to radians
    Omega_g = np.deg2rad(Omega_g)
    Fi_gm = np.deg2rad(Fi_gm)
    Lat = np.deg2rad(Lat)
    
    # computes the declination of the sun (d)
    d = 0.006918 - 0.399912 * np.cos(G) + 0.070247 * np.sin(G) - \
        0.006758 * np.cos(2 * G) + 0.000907 * np.sin(2 * G) - \
        0.002697 * np.cos(3 * G) + 0.00148 * np.sin(3 * G)
        
    # Equation of time (Et)
    # The MATLAB code uses a different simplified formula than standard models.
    Et = 0.017 + 0.4281 * np.cos(G) - 7.351 * np.sin(G) - 3.349 * np.cos(2 * G) - 9.731 * np.sin(2 * G)
    
    # computes the time of the day when the sun reaches its highest angle (tm)
    tm = 12 + (4 * (-Long) - Et) / 60
    
    # computes the hour angle of the sun (Omega_s)
    # Omega_s = (t - tm) / 12 * pi
    Omega_s = (t - tm) / 12 * np.pi
    
    # computes the zenith angle (Fi_s) (Equation 3.28 in De Bruin, simplified)
    # Fi_s = acos(sin(d) * sin(Lat) + cos(d) * cos(Lat) .* cos(Omega_s))
    Fi_s = np.arccos(np.sin(d) * np.sin(Lat) + np.cos(d) * np.cos(Lat) * np.cos(Omega_s))
    
    # computes the projected slope of the surface (Fi_g)
    # Fi_g = atan(tan(Fi_gm) .* cos(Omega_s - Omega_g))
    Fi_g = np.arctan(np.tan(Fi_gm) * np.cos(Omega_s - Omega_g))
    
    # computes the angle of the sun with the vector perpendicular to the surface (Fi_gs)
    # Fi_gs = Fi_s + Fi_g
    Fi_gs = Fi_s + Fi_g
    
    return Fi_s, Fi_gs, Fi_g, Omega_s