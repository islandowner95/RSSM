def get_bdy_from_parent_domain(parent_file,sub_file,itBdy):
    #parent_file = "ssm_out_d01.nc"
    #sub_file = "ssm_out_d02.nc"
    #itBdy = 0
    ### read latitude, longtude and variables from parent domain
    f_parent = nc4.Dataset(parent_file,'r+')
    grid_C_Lon_Deg_P = f_parent.variables["grid_C_Lon_Deg"][:,:] # parent
    grid_C_Lat_Deg_P = f_parent.variables["grid_C_Lat_Deg"][:,:]
    grid_U_Lon_Deg_P = f_parent.variables["grid_U_Lon_Deg"][:,:]
    grid_U_Lat_Deg_P = f_parent.variables["grid_U_Lat_Deg"][:,:]
    grid_V_Lon_Deg_P = f_parent.variables["grid_V_Lon_Deg"][:,:]
    grid_V_Lat_Deg_P = f_parent.variables["grid_V_Lat_Deg"][:,:]
    maskC_P = f_parent.variables["maskC"][:,:]
    maskU_P = f_parent.variables["maskU"][:,:]
    maskV_P = f_parent.variables["maskV"][:,:]
    nyP, nxP = np.shape(grid_C_Lon_Deg_P)
    nyU, nxU = np.shape(grid_U_Lon_Deg_P)
    nyV, nxV = np.shape(grid_V_Lon_Deg_P)
    eta_P = f_parent.variables["eta"][itBdy,:,:]
    u_P   = f_parent.variables["u"][itBdy,:,:]
    v_P   = f_parent.variables["v"][itBdy,:,:]
    eta_P[maskC_P==0] = 0
    u_P[maskU_P==0]   = 0
    v_P[maskV_P==0]   = 0
    
    ### read latitude, longtude from sub domain
    f_sub = nc4.Dataset(sub_file,'r+')
    grid_C_Lon_Deg_S = f_sub.variables["grid_C_Lon_Deg"][:,:] # sub
    grid_C_Lat_Deg_S = f_sub.variables["grid_C_Lat_Deg"][:,:]
    grid_U_Lon_Deg_S = f_sub.variables["grid_U_Lon_Deg"][:,:]
    grid_U_Lat_Deg_S = f_sub.variables["grid_U_Lat_Deg"][:,:]
    grid_V_Lon_Deg_S = f_sub.variables["grid_V_Lon_Deg"][:,:]
    grid_V_Lat_Deg_S = f_sub.variables["grid_V_Lat_Deg"][:,:]
    # get the north, south, east and weat boundary (lat, lon) by slices
    grid_C_Lon_North = grid_C_Lon_Deg_S[-1,:]
    grid_C_Lat_North = grid_C_Lat_Deg_S[-1,:]
    grid_C_Lon_South = grid_C_Lon_Deg_S[0,:]
    grid_C_Lat_South = grid_C_Lat_Deg_S[0,:]
    grid_C_Lon_East = grid_C_Lon_Deg_S[:,-1]
    grid_C_Lat_East = grid_C_Lat_Deg_S[:,-1]
    grid_C_Lon_West = grid_C_Lon_Deg_S[:,0]
    grid_C_Lat_West = grid_C_Lat_Deg_S[:,0]

    grid_U_Lon_North = grid_U_Lon_Deg_S[-1,:]
    grid_U_Lat_North = grid_U_Lat_Deg_S[-1,:]
    grid_U_Lon_South = grid_U_Lon_Deg_S[0,:]
    grid_U_Lat_South = grid_U_Lat_Deg_S[0,:]
    grid_U_Lon_East = grid_U_Lon_Deg_S[:,-1]
    grid_U_Lat_East = grid_U_Lat_Deg_S[:,-1]
    grid_U_Lon_West = grid_U_Lon_Deg_S[:,0]
    grid_U_Lat_West = grid_U_Lat_Deg_S[:,0]

    grid_V_Lon_North = grid_V_Lon_Deg_S[-1,:]
    grid_V_Lat_North = grid_V_Lat_Deg_S[-1,:]
    grid_V_Lon_South = grid_V_Lon_Deg_S[0,:]
    grid_V_Lat_South = grid_V_Lat_Deg_S[0,:]
    grid_V_Lon_East = grid_V_Lon_Deg_S[:,-1]
    grid_V_Lat_East = grid_V_Lat_Deg_S[:,-1]
    grid_V_Lon_West = grid_V_Lon_Deg_S[:,0]
    grid_V_Lat_West = grid_V_Lat_Deg_S[:,0]
 
    ### reshape the parent domain lat, lon and variables for griddata function
    lon_C = np.reshape(grid_C_Lon_Deg_P,(nyP*nxP,))
    lat_C = np.reshape(grid_C_Lat_Deg_P,(nyP*nxP,))
    eta_C = np.reshape(eta_P,(nyP*nxP,))

    lon_U = np.reshape(grid_U_Lon_Deg_P,(nyU*nxU,))
    lat_U = np.reshape(grid_U_Lat_Deg_P,(nyU*nxU,))
    u_U = np.reshape(u_P,(nyU*nxU,))

    lon_V = np.reshape(grid_V_Lon_Deg_P,(nyV*nxV,))
    lat_V = np.reshape(grid_V_Lat_Deg_P,(nyV*nxV,))
    v_V = np.reshape(v_P,(nyV*nxV,))
    
    ### interplate to boundary
    eta_N = griddata((lon_C,lat_C),eta_C,(grid_C_Lon_North,grid_C_Lat_North), method='linear')
    eta_S = griddata((lon_C,lat_C),eta_C,(grid_C_Lon_South,grid_C_Lat_South), method='linear')
    eta_E = griddata((lon_C,lat_C),eta_C,(grid_C_Lon_East,grid_C_Lat_East), method='linear')
    eta_W = griddata((lon_C,lat_C),eta_C,(grid_C_Lon_West,grid_C_Lat_West), method='linear')

    u_N = griddata((lon_U,lat_U),u_U,(grid_U_Lon_North,grid_U_Lat_North), method='linear')
    u_S = griddata((lon_U,lat_U),u_U,(grid_U_Lon_South,grid_U_Lat_South), method='linear')
    u_E = griddata((lon_U,lat_U),u_U,(grid_U_Lon_East,grid_U_Lat_East), method='linear')
    u_W = griddata((lon_U,lat_U),u_U,(grid_U_Lon_West,grid_U_Lat_West), method='linear')

    v_N = griddata((lon_V,lat_V),v_V,(grid_V_Lon_North,grid_V_Lat_North), method='linear')
    v_S = griddata((lon_V,lat_V),v_V,(grid_V_Lon_South,grid_V_Lat_South), method='linear')
    v_E = griddata((lon_V,lat_V),v_V,(grid_V_Lon_East,grid_V_Lat_East), method='linear')
    v_W = griddata((lon_V,lat_V),v_V,(grid_V_Lon_West,grid_V_Lat_West), method='linear')
    dictBdy = {}
    dictBdy["eta_N"] = eta_N
    dictBdy["eta_S"] = eta_S
    dictBdy["eta_E"] = eta_E
    dictBdy["eta_W"] = eta_W
    dictBdy["u_N"] = u_N
    dictBdy["u_S"] = u_S
    dictBdy["u_E"] = u_E
    dictBdy["u_W"] = u_W
    dictBdy["v_N"] = v_N
    dictBdy["v_S"] = v_S
    dictBdy["v_E"] = v_E
    dictBdy["v_W"] = v_W
    return dictBdy

#eta_N,eta_S,eta_E,eta_W,u_N,u_S,u_E,u_W,v_N,v_S,v_E,v_W = get_bdy_from_parent_domain(parent_file,sub_file,itBdy)
