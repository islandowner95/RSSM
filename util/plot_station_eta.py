import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import tools

latSite = 31.341735000000
lonSite = 121.575749000000


# read data from ERA-Interim
fileModel = nc4.Dataset("ssm_out_d02.nc" ,mode='r')
eta   = fileModel.variables["eta"][:,:,:]  # dimension [time,lat,lon]
eta   = eta*100 # m -> cm
grid_C_Lat_Deg = fileModel.variables["grid_C_Lat_Deg"][:,:]
grid_C_Lon_Deg = fileModel.variables["grid_C_Lon_Deg"][:,:]
maskC = fileModel.variables["maskC"][:,:]  
fileModel.close()

# find the nearest ponit of station
nt,ny,nx = np.shape(eta)
distMin = 200E3
for iy in range(ny):
    for ix in range(nx):
        if maskC[iy,ix] == 1:
            lon1 = grid_C_Lon_Deg[iy,ix]
            lat1 = grid_C_Lat_Deg[iy,ix]
            dist = tools.SphereDistance(lon1, lat1, lonSite, latSite, option="deg")
            if dist<distMin:
                distMin = dist
                ixSite = ix
                iySite = iy

eta_site = []
for it in range(nt):
    eta_site.append(eta[it,iySite,ixSite])
eta_site = np.array(eta_site)
x_data = np.arange(0,nt)
fig, ax = plt.subplots(1,1)
plt.plot(x_data,eta_site,color='blue',linewidth=2.0)
#ax.set_yticks(yticks)
#ax.set_xticks(xticks)
#ax.set_yticklabels(yticklabels)
#x.set_xticklabels(xticklabels)
ax.set_xlabel("Sea Level Elevation/cm")
ax.set_ylabel("Forcast Time(hours)")
ax.set_title("Gao Qiao")
figName = "fig_site/Gao_Qiao_0012.png"
fig.savefig(figName)
plt.close()




