def define_temp_response_biochem():
    """
    Translation of src/IO/define_temp_response_biochem.m
    Defines parameters for temperature response functions of biochemical processes.
    """
    TDP = {}

    # Vcmax parameters
    TDP['delHaV'] = 65330
    TDP['delSV'] = 485
    TDP['delHdV'] = 149250

    # Jmax parameters (not used in biochemical.m, but defined here)
    TDP['delHaJ'] = 43540
    TDP['delSJ'] = 495
    TDP['delHdJ'] = 152040

    # TPU parameters (commented out in biochemical.m)
    TDP['delHaP'] = 53100
    TDP['delSP'] = 490
    TDP['delHdP'] = 150650

    # Rd parameters
    TDP['delHaR'] = 46390
    TDP['delSR'] = 490
    TDP['delHdR'] = 150650

    # Kc, Ko, Gamma* parameters
    TDP['delHaKc'] = 79430
    TDP['delHaKo'] = 36380
    TDP['delHaT'] = 37830

    # C4/Dutta et al. Q10 parameters (used in C4 temperature correction)
    TDP['Q10'] = 2
    TDP['s1'] = 0.3
    TDP['s2'] = 313.15
    TDP['s3'] = 0.2
    TDP['s4'] = 288.15
    TDP['s5'] = 1.3
    TDP['s6'] = 328.15
    
    return TDP