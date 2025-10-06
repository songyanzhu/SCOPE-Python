import numpy as np

def Soil_Inertia0(cs, rhos, lambdas):
    """
    Calculates soil thermal inertia using fundamental soil properties.
    
    :param cs: Soil specific heat capacity [J kg-1 K-1].
    :param rhos: Soil bulk density [kg m-3].
    :param lambdas: Soil thermal conductivity [W m-1 K-1].
    :return: Soil thermal inertia GAM [J m-2 K-1 s-0.5].
    """
    # GAM = sqrt(cs * rhos * lambdas)
    GAM = np.sqrt(cs * rhos * lambdas)
    return GAM