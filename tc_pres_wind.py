
def get_ty_pres(gridX,gridY,Lon0,Lat0,Xc,Yc,Po,Rm):
    """
    Po: pressure at TC center, Pa
    Xc, Yc: TC center position
    output: return "windU" or "windV"
    """
    rCX     = tools.SphereDistance(Lon0,Yc,Xc,Yc,option="deg")  # TC distance from [0,0]
    rCY     = tools.SphereDistance(Xc,Lat0,Xc,Yc,option="deg")
    r       = np.zeros(np.shape(gridX)) # distance from TC center
    r       = np.sqrt( (gridX-rCX)**2 + (gridY-rCY)**2 )
    r[r==0] = 1.0 # avoild division by zero
    deltaP  = Pe-Po
    if ty_pres_wind == "Fujita-Takah":
        Pres    = Pe - deltaP/np.sqrt(1+2*(r/Rm)**2) # Fujita     0<r<2Rm
        Pres2   = Pe - deltaP/(1+r/Rm)                 # Takahashi  2Rm<=
        Pres[r>=2*Rm] = Pres2[r>=2*Rm]
    elif ty_pres_wind == "Fujita":
        Pres    = Pe - deltaP/np.sqrt(1+2*(r/Rm)**2) # Fujita  
    elif ty_pres_wind == "Takah":
        Pres   = Pe - deltaP/(1+r/Rm)                 # Takahashi  
    elif ty_pres_wind == "Holland":
        B = 1.5 + ((98000-101000+deltaP)/100)/120
        #B       = 1.05 # Holland B parameter
        Pres    = Po + deltaP*np.exp(-1.0*(Rm/r)**B) # Holland
    else:
        print("Error: the ty_pres_wind option is not exsit")
        exit()
    return Pres

def get_ty_wind_stress(gridX,gridY,grid_Lat,Lon0,Lat0,Vcx,Vcy,Xc,Yc,Po,Rm,output):
    """
    gridX, gridY: Cartesian Coord X and Y, meter
    Lon0, Lat0: rdinate origin location, degree
    Po: pressure at TC center, Pa
    Vcx, Vcy: TC transfer velocity, m/s
    Xc, Yc: TC center position
    output: return "windU" or "windV"
    """
    inflowA = 20/180*math.pi # inflow angle, radian
    rCX     = tools.SphereDistance(Lon0,Yc,Xc,Yc,option="deg")  # TC distance from [0,0]
    rCY     = tools.SphereDistance(Xc,Lat0,Xc,Yc,option="deg")
    r       = np.zeros(np.shape(gridX)) # distance from TC center
    r       = np.sqrt( (gridX-rCX)**2 + (gridY-rCY)**2 )
    distPX  = gridY - rCY # distance betwwn grid Position and X-axis of TC center
    distPY  = gridX - rCX # distance betwwn grid Position and Y-axis of TC center
    theta0  = np.arcsin(distPX/r) # angle between CP and CX, C is the center of typhoon, P is the points of grid, CX is the X-axis, radian, 0-3.1415926
    theta0[distPY<0] = math.pi-theta0[distPY<0]
    theta0[theta0<0] = 2*math.pi + theta0[theta0<0]
    theta   = theta0 + inflowA
    deltaP = Pe-Po
    f = 2.0*omega*np.sin(grid_Lat)
    r[r==0] = 1.0 # avoild division by zero
    # gradient wind
    if ty_pres_wind == "Fujita-Takah":
        Vgrad  = np.sqrt( 2*deltaP/rho_a*(r/Rm)**2*( (1+2*(r/Rm)**2)**-1.5 ) + (r*f/2.0)**2 ) - r*f/2.0 # Fujita 0<r<2R
        Vgrad2 = np.sqrt( deltaP/rho_a*r/Rm/(1+r/Rm)**2 + (r*f/2.0)**2 ) - r*f/2.0                      # Takahashi  2Rm<=
        Vgrad[r>=2*Rm] = Vgrad2[r>=2*Rm]
    elif ty_pres_wind == "Fujita":
        Vgrad = np.sqrt( 2*deltaP/rho_a*(r/Rm)**2*( (1+2*(r/Rm)**2)**-1.5 ) + (r*f/2.0)**2 ) - r*f/2.0
    elif ty_pres_wind == "Takah":
        Vgrad = np.sqrt( deltaP/rho_a*r/Rm/(1+r/Rm)**2 + (r*f/2.0)**2 ) - r*f/2.0
    elif ty_pres_wind == "Holland":
        #B = 1.05 # Holland B parameter # wind speed intend to smaller
        B = 1.5 + ((98000-101000+deltaP)/100)/120
        Vgrad = np.sqrt( deltaP*B/rho_a*(Rm/r)**B*np.exp(-1.0*(Rm/r)**B) + (r*f/2.0)**2 ) - r*f/2.0# gradient wind
    else:
        print("Error: the ty_pres_wind option is not exsit")
        exit()
    Vgx   = -1.0*Vgrad*np.sin(theta) # x-component of Vgrad, m/s
    Vgy   = Vgrad*np.cos(theta)  # y-componest of Vgrad, m/s
    # typhoon basic wind field, Ueno Takeo(1981)
    arg0 = np.exp(-np.pi/4.0*np.abs(r-Rm)/Rm)
    Vbx = arg0*Vcx
    Vby = arg0*Vcy
    # compose wind( TC transfer velocity + gradient wind)
    c1 = 1.0 ; c2 = 0.8
    Utc = c1*Vbx + c2*Vgx
    Vtc = c1*Vby + c2*Vgy
    # wind stress
    UVtc = np.sqrt(Utc**2+Vtc**2) # total wind, m/s
    Cd = (0.8+0.065*UVtc)/1000 # drag coefficient
    if output == 'windU':
        wind_stress = rho_a*Cd*Utc*UVtc
    elif output == 'windV':
        wind_stress = rho_a*Cd*Vtc*UVtc
    else:
        print("ERROR: output option error")
        exit()
    return wind_stress


