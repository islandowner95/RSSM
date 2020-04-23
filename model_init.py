import numpy as np
from namelist_python import read_namelist_file
import tools
from scipy.interpolate import griddata
import netCDF4 as nc4
import shutil
import os
from skimage import measure

## import constants
path = os.getcwd() + '/'
exec(open(path+'constants.py').read())
get_constants()
exec(open(path+'read_namelist.py').read())
read_namelist()


def generate_grid(lonMin,lonMax,latMin,latMax,dx,dy,output_file): 
    """
    input: lonMax,lonMin,latMax,latMin,dx,dy unit is degree
    output: Simulation grid,numpy array(ny,nx)
    """
    ### grid setting
    ny = int((latMax-latMin)/dx)
    nx = int((lonMax-lonMin)/dy)
    print("nx, ny = ",nx, ny)
    Nx_c = nx # grid size
    Ny_c = ny
    Nx_u = nx+1
    Ny_u = ny
    Nx_v = nx
    Ny_v = ny+1
    
    grid_C_Lat_Deg = np.zeros((Ny_c,Nx_c))
    grid_C_Lon_Deg = np.zeros((Ny_c,Nx_c))
    grid_U_Lat_Deg = np.zeros((Ny_u,Nx_u))
    grid_U_Lon_Deg = np.zeros((Ny_u,Nx_u))
    grid_V_Lat_Deg = np.zeros((Ny_v,Nx_v))
    grid_V_Lon_Deg = np.zeros((Ny_v,Nx_v))

    for j in range(ny):
        grid_C_Lat_Deg[j,:] = latMin + dy*j
        grid_U_Lat_Deg[j,:] = latMin + dy*j
    for i in range(nx):
        grid_C_Lon_Deg[:,i] = lonMin + dx*i
        grid_V_Lon_Deg[:,i] = lonMin + dx*i
    for i in range(nx+1):
        grid_U_Lon_Deg[:,i] = lonMin - dx/2.0 + dx*i
    for j in range(ny+1):
        grid_V_Lat_Deg[j,:] = latMin - dy/2.0 + dy*j
    
    # degree to radian
    grid_C_Lat = grid_C_Lat_Deg * deg2rad
    grid_C_Lon = grid_C_Lon_Deg * deg2rad
    grid_U_Lat = grid_U_Lat_Deg * deg2rad
    grid_U_Lon = grid_U_Lon_Deg * deg2rad
    grid_V_Lat = grid_V_Lat_Deg * deg2rad
    grid_V_Lon = grid_V_Lon_Deg * deg2rad
    return_dict = dict()
    return_dict["grid_C_Lon"] = grid_C_Lon
    return_dict["grid_C_Lat"] = grid_C_Lat
    return_dict["grid_U_Lon"] = grid_U_Lon
    return_dict["grid_U_Lat"] = grid_U_Lat
    return_dict["grid_V_Lon"] = grid_V_Lon
    return_dict["grid_V_Lat"] = grid_V_Lat
    return_dict["grid_C_Lon_Deg"] = grid_C_Lon_Deg
    return_dict["grid_C_Lat_Deg"] = grid_C_Lat_Deg
    return_dict["grid_U_Lon_Deg"] = grid_U_Lon_Deg
    return_dict["grid_U_Lon_Deg"] = grid_U_Lon_Deg
    return_dict["grid_V_Lat_Deg"] = grid_V_Lat_Deg
    return_dict["grid_V_Lat_Deg"] = grid_V_Lat_Deg
  
    f = nc4.Dataset(output_file,'w',format='NETCDF4')
    f.createDimension('lon',len(grid_C_Lon[0,:])) 
    f.createDimension('lat',len(grid_C_Lat[:,0]))
    f.createDimension('lonStag',len(grid_U_Lon[0,:])) 
    f.createDimension('latStag',len(grid_V_Lat[:,0]))
    lonC = f.createVariable("lon",'f4',("lon")) 
    latC = f.createVariable("lat",'f4',("lat")) 
    lonU = f.createVariable("lonStag",'f4',("lonStag")) 
    latV = f.createVariable("latStag",'f4',("latStag")) 
    gridCLon = f.createVariable('grid_C_Lon','f8',('lat','lon'))
    gridCLat = f.createVariable('grid_C_Lat','f8',('lat','lon'))
    gridULon = f.createVariable('grid_U_Lon','f8',('lat','lonStag'))
    gridULat = f.createVariable('grid_U_Lat','f8',('lat','lonStag'))
    gridVLon = f.createVariable('grid_V_Lon','f8',('latStag','lon'))
    gridVLat = f.createVariable('grid_V_Lat','f8',('latStag','lon'))
    gridCLonDeg = f.createVariable('grid_C_Lon_Deg','f8',('lat','lon'))
    gridCLatDeg = f.createVariable('grid_C_Lat_Deg','f8',('lat','lon'))
    gridULonDeg = f.createVariable('grid_U_Lon_Deg','f8',('lat','lonStag'))
    gridULatDeg = f.createVariable('grid_U_Lat_Deg','f8',('lat','lonStag'))
    gridVLonDeg = f.createVariable('grid_V_Lon_Deg','f8',('latStag','lon'))
    gridVLatDeg = f.createVariable('grid_V_Lat_Deg','f8',('latStag','lon'))
    lonC.units = 'degrees_east'
    lonC.long_name = 'longitude'
    latC.units = 'degrees_north'
    latC.long_name = 'latitude' 
    lonU.units = 'degrees_east'
    lonU.long_name = 'longitude'
    latV.units = 'degrees_north'
    latV.long_name = 'latitude' 

    gridCLon.units = 'radian_east'
    gridCLon.long_name = 'longitude'
    gridCLat.units = 'radian_north'
    gridCLat.long_name = 'latitude' 
    gridULon.units = 'radian_east'
    gridULon.long_name = 'longitude'
    gridULat.units = 'radian_north'
    gridULat.long_name = 'latitude' 
    gridVLon.units = 'radian_east'
    gridVLon.long_name = 'longitude'
    gridVLat.units = 'radian_north'
    gridVLat.long_name = 'latitude' 

    gridCLonDeg.units = 'degrees_east'
    gridCLonDeg.long_name = 'longitude'
    gridCLatDeg.units = 'degrees_north'
    gridCLatDeg.long_name = 'latitude' 
    gridULonDeg.units = 'degrees_east'
    gridULonDeg.long_name = 'longitude'
    gridULatDeg.units = 'degrees_north'
    gridULatDeg.long_name = 'latitude' 
    gridVLonDeg.units = 'degrees_east'
    gridVLonDeg.long_name = 'longitude'
    gridVLatDeg.units = 'degrees_north'
    gridVLatDeg.long_name = 'latitude' 
   
    #depthV._FillValue = -1.0e+20
    lonC[:]   = grid_C_Lon_Deg[0,:] 
    latC[:]   = grid_C_Lat_Deg[:,0] 
    lonU[:]   = grid_U_Lon_Deg[0,:] 
    latV[:]   = grid_V_Lat_Deg[:,0] 
    gridCLon[:,:] = grid_C_Lon
    gridCLat[:,:] = grid_C_Lat
    gridULon[:,:] = grid_U_Lon
    gridULat[:,:] = grid_U_Lat
    gridVLon[:,:] = grid_V_Lon
    gridVLat[:,:] = grid_V_Lat
    gridCLonDeg[:,:] = grid_C_Lon_Deg
    gridCLatDeg[:,:] = grid_C_Lat_Deg
    gridULonDeg[:,:] = grid_U_Lon_Deg
    gridULatDeg[:,:] = grid_U_Lat_Deg
    gridVLonDeg[:,:] = grid_V_Lon_Deg
    gridVLatDeg[:,:] = grid_V_Lat_Deg
    f.close()

    return return_dict 


def inject_bathymetry(grid_Lon_Deg,grid_Lat_Deg,etopo_file,grid_file):
    """
    input:
         grid: latitude and longitude, unit: degree
         etopo_file: file name of ETOPO, NETCDF
         output_file: output file name, NETCDF
    output: none
    """
    etopo = nc4.Dataset(etopo_file,'r+')
    lat0  = etopo.variables['lat']
    lon0  = etopo.variables['lon']
    topo0 = etopo.variables['Band1']
    nx = len(lon0); ny = len(lat0)
    lon1,lat1 = np.meshgrid(lon0,lat0)
    lon2  = np.reshape(lon1,(nx*ny,)) 
    lat2  = np.reshape(lat1,(nx*ny,)) 
    topo2 = np.reshape(topo0,(nx*ny,)) 
    
    interpMethod = "nearest"
    #interpMethod = "linear"
    print("interpolating bathymetry")
    depthC = griddata((lon2,lat2), topo2,(grid_Lon_Deg,grid_Lat_Deg), method=interpMethod)
    depthC = -1.0*depthC # depth below sea level transfer as positive, depth > 0 as water, depth < 0 as land
    """
    depthC = np.zeros(np.shape(grid_Lon_Deg))
    Nyy,Nxx = np.shape(depthC)
    for i in range(Nxx):
        depthC[:,i] = 5+i/Nxx*15
    depthC[[0,-1],:] = 0
    depthC[:,[0,-1]] = 0
    """
    f = nc4.Dataset(grid_file,'a',format='NETCDF4')
    depth = f.createVariable('depth','f8',('lat','lon'))
    depth.units = 'm'
    depth.long_name = 'bathymetry' 
    #depthV._FillValue = -1.0e+20
    depth[:,:] = depthC  
    f.close()
    return depthC

def mark_land_bdy(grid_file):
    """
    ocean mask = 1
    land and inland lake mask = 0
    """
    print("Marking the ocean and land cell...")
    f1 = nc4.Dataset(grid_file,'r')
    depthC = f1.variables["depth"][:,:]
    lonC   = f1.variables["lon"][:]
    latC   = f1.variables["lat"][:]
    f1.close()
    
    maskC  = np.zeros(np.shape(depthC))
    nx = len(lonC)
    ny = len(latC)
    maskC[depthC>5.0] = 1    # depth > 5m as water, ocean and lake = 1, land = 0
    maskC = cull_lake(maskC) # ocean = 1, land and inland lake = 0

    maskU = np.ones((ny,nx+1)) # u mask
    maskV = np.ones((ny+1,nx)) # v mask
    for j in range(0,ny):
        for i in range(1,nx):
            maskU[j,i] = maskC[j,i-1]*maskC[j,i]
    for j in range(1,ny):
        for i in range(0,nx):
            maskV[j,i] = maskC[j-1,i]*maskC[j,i]
    maskU[:,[0,-1]] = maskC[:,[0,-1]]
    maskV[[0,-1],:] = maskC[[0,-1],:] 

    f2 = nc4.Dataset(grid_file,mode='a')
    maskC1 = f2.createVariable('maskC','f8',('lat','lon'))
    maskU1 = f2.createVariable('maskU','f8',('lat','lonStag'))
    maskV1 = f2.createVariable('maskV','f8',('latStag','lon'))
    maskC1.long_name = 'mask at Center position'
    maskU1.long_name = 'mask at U position'
    maskV1.long_name = 'mask at V position'
    maskC1[:,:] = maskC 
    maskU1[:,:] = maskU 
    maskV1[:,:] = maskV 
    f2.close()
    return None   

def cull_lake(maskOld):
    """
    input:maskOld include lake
    output:maskNew exclude lake
    maskOld: value = 0, depth>0, land, as background
             value = 1, depth<0, ocean or lake

    maskNew: value = 0, land or inland lake
             value = 1, ocean
    """
    #area_labels = measure.label(maskOld, connectivity=2)
    area_labels = measure.label(maskOld, connectivity=1)
    target_index = 0 # ocean
    count_max = 0
    max_index = np.max(area_labels) 
    for index in range(1,max_index+1):
        count = np.sum(area_labels==index)
        if count > count_max:
            count_max = count
            target_index = index
    #print(target_index,count_max)
    maskNew = np.zeros(np.shape(maskOld))
    maskNew[area_labels!=target_index] = 2 # mark land and inland lake = 2
    maskNew[maskNew==0] = 1 # mark land and inland lake = 1
    maskNew[maskNew==2] = 0 # mark land and inland lake = 0
    return maskNew

def get_tc_best_track():
    """0
    get typhoon best track
    Typhoon: lat, lon, Rmax, P0
    """
    print("Proessing typhoon best track data...") 
    namelist   = read_namelist_file('namelist.ssm')
    # namelist 
    startDate = namelist.groups["share"]["startDate"]
    endDate   = namelist.groups["share"]["endDate"]
    interval  = namelist.groups["typhoon"]["ty_interval"]
    tcBSTFile = namelist.groups["typhoon"]["tyBSTFile"]
    f = nc4.Dataset(tcBSTFile,mode='w')
    # creat time dimension
    timeNum = tools.get_time_dimension(startDate,endDate,interval)
    f.createDimension('times', timeNum)    
    # creat variables
    time  = f.createVariable('times','S19',('times'))
    Po    = f.createVariable('Po','f8',('times'))
    RMW   = f.createVariable('RMW','f8',('times'))
    lonTy = f.createVariable('lonTy','f8',('times'))
    latTy = f.createVariable('latTy','f8',('times'))
    VT    = f.createVariable('VT','f8',('times'))
    Theta = f.createVariable('Theta','f8',('times'))
    time.units = 'yyyy-mm-dd_hh:ss:mm'
    time.long_name = 'time of typhoon records'
    time.calendar = "gregorian"
    Po.units = 'Pa'
    Po.long_name = 'Pressure at typhoon center'
    RMW.units = 'm'
    RMW.long_name = 'Typhoon maximum wind radius'
    lonTy.units = 'degrees_east'
    lonTy.long_name = 'Longitude of Typhoon center'
    latTy.units = 'degrees_north'
    latTy.long_name = 'Latitude of Typhoon center'
    VT.units = 'm/s'
    VT.long_name = 'typhoon transfer volocity'
    Theta.units = 'degree'
    Theta.long_name = 'typhoon transfer direction'
    dictTy = tools.read_typhoon_record()
    print(dictTy["dateTy"])
    time  = dictTy["dateTy"]
    Po[:]    = dictTy["presTy"]
    RMW[:]   = dictTy["rmwTy"]
    lonTy[:] = dictTy["lonTy"]
    latTy[:] = dictTy["latTy"]
    VT[:]    = dictTy["VT"]
    Theta[:] = dictTy["Theta"]
    f.close()
    print("Output file:",tcBSTFile)
    return None    

def gen_init_bdy_from_model():
    """
    Generating presure gradiant force and wind stress from ERA-Interim data.
    """
    namelist  = read_namelist_file('namelist.ssm')
    dx = namelist.groups["domain1"]["dx"]
    dxRad = tools.deg2rad(dx)
    R = namelist.groups["constant"]["R"]*1000 # km -> m
    # reading ERA-Interim data
    fileModel = nc4.Dataset("ERA-Interim201907.nc" ,mode='r')
    u10 = fileModel.variables["u10"][:,:,:] #[time,lat,lon]
    #u10Scale = fileModel.variables["u10"].scale_factor
    #u10Add   = fileModel.variables["u10"].add_offset
    #u10 = u10*u10Scale+u10Add
    v10 = fileModel.variables["u10"][:,:,:] #[time,lat,lon]
    msl = fileModel.variables["u10"][:,:,:] #[time,lat,lon]
    lat = fileModel.variables["latitude"][:]  #[lat]
    lon = fileModel.variables["longitude"][:] #[lon]
    nx = len(lon); ny = len(lat)
    lon1,lat1 = np.meshgrid(lon,lat)
    lon2  = np.reshape(lon1,(nx*ny,))
    lat2  = np.reshape(lat1,(nx*ny,))
    fileModel.close()
    # domain grid information
    initBdyD01 = namelist.groups["domain1"]["initBdyD01"]
    fileInitBdy = nc4.Dataset(initBdyD01,mode='a')
    lonC = fileInitBdy.variables["lonC"][:] 
    latC = fileInitBdy.variables["latC"][:] 
    lonU = fileInitBdy.variables["lonU"][:] 
    latU = fileInitBdy.variables["latU"][:] 
    lonV = fileInitBdy.variables["lonV"][:] 
    latV = fileInitBdy.variables["latV"][:] 
    nxC = len(lonC)
    nyC = len(latC)
    nxU = len(lonU)
    nyU = len(latU)
    nxV = len(lonV)
    nyV = len(latV)
    gridXC,gridYC = np.meshgrid(lonC, latC)
    gridXU,gridYU = np.meshgrid(lonU, latU)
    gridXV,gridYV = np.meshgrid(lonV, latV)
    gridYURad = tools.deg2rad(gridYU)
    gridYVRad = tools.deg2rad(gridYV)
    lonU2 = np.zeros([nxU+2])
    lonU2[1:-1] = lonU
    lonU2[0]    = lonU[0]-0.5*(lonU[1]-lon[0])
    lonU2[-1]   = lonU[-1]+0.5*(lonU[-1]-lonU[-2])
    gridXU2,gridYU2 = np.meshgrid(lonU2, latU)
    latV2 = np.zeros([nyV+2])
    latV2[1:-1] = latV
    latV2[0]    = latV[0]-0.5*(latV[1]-latV[0])
    latV2[-1]   = latV[-1]+0.5*(latV[-1]-latV[-2])
    gridXV2,gridYV2 = np.meshgrid(lonV, latV2)
    # interplate
    #time = range(68,77) 
    #time = range(68,77) 
    time = range(60,69) 
    nt = len(time)
    presC = np.zeros([nt,nyC,nxC]) # pressure at cell center position / Pa
    presU = np.zeros([nyU,nxU+2])  # pressure at U position to calculate PGFu
    presV = np.zeros([nyV+2,nxV])  # pressure at V position to calculate PGFv
    u10U = np.zeros([nt,nyU,nxU]) # u10 and v10 at U position
    v10U = np.zeros([nt,nyU,nxU])
    Fs   = np.zeros([nt,nyU,nxU]) # x direction wind stress at U position
    PGFu = np.zeros([nt,nyU,nxU]) # x direction pressure gradiant force  at U position
    u10V = np.zeros([nt,nyV,nxV]) # u10 and v10 at U position
    v10V = np.zeros([nt,nyV,nxV])
    Gs   = np.zeros([nt,nyV,nxV]) # y direction wind stress at V position
    PGFv = np.zeros([nt,nyV,nxV]) # y direction pressure gradiant force  at V position
    for i in range(len(time)): 
        it = time[i] # ERA data time dimension 
        ### wind stress    
        u10_i = np.reshape(u10[it,:,:],(nx*ny,))
        v10_i = np.reshape(v10[it,:,:],(nx*ny,))
        u10U[i,:,:]  = griddata((lon2,lat2), u10_i,(gridXU,gridYU), method='linear')
        v10U[i,:,:]  = griddata((lon2,lat2), v10_i,(gridXU,gridYU), method='linear')
        u10V[i,:,:]  = griddata((lon2,lat2), u10_i,(gridXV,gridYV), method='linear')
        v10V[i,:,:]  = griddata((lon2,lat2), v10_i,(gridXV,gridYV), method='linear')
        Fs[i,:,:]    = para.u_wind_stress(u10U[i,:,:],v10U[i,:,:])  # x-direction wind stress at U position
        Gs[i,:,:]    = para.v_wind_stress(u10V[i,:,:],v10V[i,:,:])  # y-direction wind stress at V position
        ### Pressure and PGF
        msl_i = np.reshape(msl[it,:,:],(nx*ny,))
        presC[i,:,:] = griddata((lon2,lat2), msl_i,(gridXC,gridYC), method='linear') 
        presU[:,:]   = griddata((lon2,lat2), msl_i,(gridXU2,gridYU2), method='linear') 
        presV[:,:]   = griddata((lon2,lat2), msl_i,(gridXV2,gridYV2), method='linear') 
        # PGFu
        argU = 1/(R*np.cos(gridYURad*2*dxRad))
        PGFu[i,:,:] = argU*(presU[:,2:]-presU[:,0:-2])
        # PGFv
        argV = 1/(R*2*dxRad)
        PGFv[i,:,:] = argV*(presV[2:,:]-presV[0:-2,:])
    # save data
    PGFu0 = fileInitBdy.createVariable('PGFu','f8',('time','lat','lonStag'))
    PGFu0.units = 'N'
    PGFu0.long_name = 'pressure gradiant force at U position' 
    PGFv0 = fileInitBdy.createVariable('PGFv','f8',('time','latStag','lon'))
    PGFv0.units = 'N'
    PGFv0.long_name = 'pressure gradiant force at V position' 
    Fs0 = fileInitBdy.createVariable('Fs','f8',('time','lat','lonStag'))
    Fs0.units = ''
    Fs0.long_name = 'x-direction wind stress at U position' 
    Gs0 = fileInitBdy.createVariable('Gs','f8',('time','latStag','lon'))
    Gs0.units = ''
    Gs0.long_name = 'y-direction wind stress at V position' 
    presC0 = fileInitBdy.createVariable('presC','f8',('time','lat','lon'))
    presC0.units = 'Pa'
    presC0.long_name = 'pressure at cell center' 
    print(np.shape(PGFu))
    print(np.shape(PGFu0))
    PGFu0[:,:,:] = PGFu
    PGFv0[:,:,:] = PGFv
    Fs0[:,:,:] = Fs
    Gs0[:,:,:] = Gs
    presC0[:,:,:] = presC
    fileInitBdy.close() 
    print("Generating presure gradiant force and wind stress from ERA-Interim data")
    return None

if __name__ == '__main__':
    domain1 = False
    domain2 = False
     
    if domain1:
        ID = dom_id[0]-1
        print("domina ",ID)
        # generate grid and inject depth
        dxDeg = dx[ID]*rad2deg
        dyDeg = dy[ID]*rad2deg
        return_dict = generate_grid(lonMin[ID],lonMax[ID],latMin[ID],latMax[ID],dxDeg,dyDeg,gridFile[ID])
        grid_C_Lon_Deg = return_dict["grid_C_Lon_Deg"]
        grid_C_Lat_Deg = return_dict["grid_C_Lat_Deg"]
        depthC = inject_bathymetry(grid_C_Lon_Deg,grid_C_Lat_Deg,etopoFile,gridFile[ID])
        print(np.min(depthC),np.max(depthC))
        # mark land and inland lake
        mark_land_bdy(gridFile[ID])

    if domain2:
        ID = dom_id[1]-1 
        print("domina ",ID)
        # generate grid and inject depth
        dxDeg = dx[ID]*rad2deg
        dyDeg = dy[ID]*rad2deg
        return_dict = generate_grid(lonMin[ID],lonMax[ID],latMin[ID],latMax[ID],dxDeg,dyDeg,gridFile[ID])
        grid_C_Lon_Deg = return_dict["grid_C_Lon_Deg"]
        grid_C_Lat_Deg = return_dict["grid_C_Lat_Deg"]
        depthC = inject_bathymetry(grid_C_Lon_Deg,grid_C_Lat_Deg,etopoFile,gridFile[ID])
        print(np.min(depthC),np.max(depthC))
        # mark land and inland lake
        mark_land_bdy(gridFile[ID])

    # get tc best track
    get_tc_best_track()
    # Generating presure gradiant force and wind stress from ERA-Interim data
    #gen_init_bdy_from_model()     








