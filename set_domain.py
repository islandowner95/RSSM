def set_domain(domain_id):
    global grid_C_Lon_Deg,grid_C_Lat_Deg
    global grid_U_Lon_Deg,grid_U_Lat_Deg
    global grid_V_Lon_Deg,grid_V_Lat_Deg
    global grid_C_Lon,grid_C_Lat
    global grid_U_Lon,grid_U_Lat
    global grid_V_Lon,grid_V_Lat
    global maskC,maskU,maskV
    global H
    global nx,ny,Nx_c,Ny_c,Nx_u,Ny_u,Nx_v,Ny_v
    global f_U,f_V
    global gridXC, gridYC
    global gridXU, gridYU
    global gridXV, gridYV
    ID = domain_id - 1
    f_grid  = nc4.Dataset(gridFile[ID],'r')
    # grid of domain1 
    grid_C_Lon_Deg = f_grid.variables["grid_C_Lon_Deg"][:,:]
    grid_C_Lat_Deg = f_grid.variables["grid_C_Lat_Deg"][:,:]
    grid_U_Lon_Deg = f_grid.variables["grid_U_Lon_Deg"][:,:]
    grid_U_Lat_Deg = f_grid.variables["grid_U_Lat_Deg"][:,:]
    grid_V_Lon_Deg = f_grid.variables["grid_V_Lon_Deg"][:,:]
    grid_V_Lat_Deg = f_grid.variables["grid_V_Lat_Deg"][:,:]
    grid_C_Lon     = f_grid.variables["grid_C_Lon"][:,:]
    grid_C_Lat     = f_grid.variables["grid_C_Lat"][:,:]
    grid_U_Lon     = f_grid.variables["grid_U_Lon"][:,:]
    grid_U_Lat     = f_grid.variables["grid_U_Lat"][:,:]
    grid_V_Lon     = f_grid.variables["grid_V_Lon"][:,:]
    grid_V_Lat     = f_grid.variables["grid_V_Lat"][:,:]
    
    # depth H and mask of domain1
    maskC = f_grid.variables["maskC"][:,:]
    maskU = f_grid.variables["maskU"][:,:]
    maskV = f_grid.variables["maskV"][:,:]
    H     = f_grid.variables["depth"][:,:]# depth of water
    H     = H*maskC
    H[H>500] = 500

    # dimension of domain1
    ny, nx = np.shape(maskC)
    Nx_c = nx # grid size
    Ny_c = ny
    Nx_u = nx+1
    Ny_u = ny
    Nx_v = nx
    Ny_v = ny+1

    ### coriolis force parameter
    f_U = 2*omega*np.sin(grid_U_Lat)
    f_V = 2*omega*np.sin(grid_V_Lat)

    # Cartesian Coord of domain1 grid
    gridXC, gridYC = getCartesianCoord(grid_C_Lon,grid_C_Lat)
    gridXU, gridYU = getCartesianCoord(grid_U_Lon,grid_U_Lat)
    gridXV, gridYV = getCartesianCoord(grid_V_Lon,grid_V_Lat)
    
    return None

def getCartesianCoord(grid_Lon,grid_Lat):
    NY, NX = np.shape(grid_Lon)
    gridX = np.zeros([NY,NX])
    gridY = np.zeros([NY,NX])
    for j in range(NY):
        for i in range(NX):
            lon1 = grid_Lon[j,i]
            lat1 = grid_Lat[j,i]

            lon0 = grid_Lon[0,0]
            lat0 = lat1
            gridX[j,i] = tools.SphereDistance(lon0,lat0,lon1,lat1,option="rad")

            lon0 = lon1
            lat1 = grid_Lat[0,0]
            gridY[j,i] = tools.SphereDistance(lon0,lat0,lon1,lat1,option="rad")
    return gridX, gridY


