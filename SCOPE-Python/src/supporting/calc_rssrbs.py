def calc_rssrbs(SMC, LAI, rbs):
    """
    Calculates soil surface resistance (rss) and bulk resistance (rbs).
    
    :param SMC: Soil Moisture Content [m3/m3].
    :param LAI: Leaf Area Index.
    :param rbs: Initial bulk resistance for the canopy [s m-1].
    :return: A tuple (rss, rbs) of calculated resistances.
    """
    
    # rss = 11.2 * exp(42 * (0.22 - SMC))
    import numpy as np
    rss = 11.2 * np.exp(42 * (0.22 - SMC))
    
    # rbs = rbs * LAI / 3.3
    # Note: This overwrites the input `rbs` value
    rbs = rbs * LAI / 3.3
    
    return rss, rbs