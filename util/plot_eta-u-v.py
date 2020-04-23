import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt


# read data from ERA-Interim
fileModel = nc4.Dataset("../ssm_out_d02.nc" ,mode='r')
u     = fileModel.variables["u"][:,:,:]  # dimension [time,lat,lon]
v     = fileModel.variables["v"][:,:,:]  # dimension [time,lat,lon]
eta   = fileModel.variables["eta"][:,:,:]  # dimension [time,lat,lon]
eta   = eta*100 # m -> cm
maskU = fileModel.variables["maskU"][:,:]  
maskV = fileModel.variables["maskV"][:,:]  
maskC = fileModel.variables["maskC"][:,:]  
grid_C_Lon_Deg = fileModel.variables["grid_C_Lon_Deg"][:,:]  
grid_C_Lat_Deg = fileModel.variables["grid_C_Lat_Deg"][:,:]  
lon = fileModel.variables["lon"][:]
lat = fileModel.variables["lat"][:]
fileModel.close()
nt,ny,nx = np.shape(u) # get the dimension
for it in range(nt):
    u[it,maskU==0]     = 0
    v[it,maskV==0]     = 0
    eta[it,maskC==0]   = 0

#eta_max = 300
#eta_min = -100
eta_max = np.max(eta)
eta_min = np.min(eta)
levels=np.linspace(-20,100,13)
for i in range(nt):
    fig, ax = plt.subplots(1,1)
    eta_p = eta[i,:,:]
    eta_p[maskC==0] = np.nan
    #plt_eta = ax.contourf(grid_C_Lon_Deg,grid_C_Lat_Deg,eta_p,cmap='jet')
    plt_eta = ax.contourf(grid_C_Lon_Deg,grid_C_Lat_Deg,eta_p,cmap='jet',levels=levels)
    cbar = plt.colorbar(plt_eta,label='cm')
    cbar.set_ticks(np.linspace(-20,100,13))
 
    xx = []; yy = []; uu = []; vv = []
    u_eta = (u[i,:,1:]+u[i,:,:-1])/2.0 # u at eta position
    v_eta = (v[i,1:,:]+v[i,:-1,:])/2.0 # v at eta position
    u_eta[maskC==0] = np.nan
    v_eta[maskC==0] = np.nan
    skip = 8
    for irow in range( 0, ny, skip ):
        for icol in range( 0, nx, skip ):
            iLon = lon[icol]
            iLat = lat[irow]
            xx.append(iLon); yy.append(iLat)
            uu.append( u_eta[irow,icol] )
            vv.append( v_eta[irow,icol] )
    quiv = ax.quiver( xx, yy, uu, vv, color='blue', scale=1.0,scale_units="inches")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title("Forecast Time: "+str(i)+" Hours")
    figName = "fig2/eta-u-v_"+str(i)+"h.png"
    fig.savefig(figName)
    plt.close()




