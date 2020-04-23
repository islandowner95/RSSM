def get_constants():
    global R,omega,rho_a,rho_w,Pe,deg2rad,rad2deg,g,dragConst,Ah
    R         =  6378.135E3  # Earth radius
    omega     = 7.292E-05    #
    rho_a     = 1.293        # Air density, kg*m-3
    rho_w     = 1025.0       # Water density, kg*m-3
    Pe        = 1010.0E2     # Ambient pressure of Typhoon, Pa
    deg2rad   = 3.1415926/180.0 # 
    rad2deg   = 180.0/3.1415926 # 
    g         = 9.81     # gravity accelation, m/s
    #dragConst = 2.6E-3
    dragConst = 1.0E-3
    Ah        = 2.5
    return None
