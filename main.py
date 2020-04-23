import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import matplotlib as mpl
import tools
import os
import shutil
import time
from namelist_python import read_namelist_file
import netCDF4 as nc4
from scipy.interpolate import griddata

start_time = time.time()

## import constants and namelist
path = os.getcwd() + '/'
exec(open(path+'constants.py').read())
get_constants()
exec(open(path+'read_namelist.py').read())
read_namelist()
exec(open(path+'set_domain.py').read())
exec(open(path+'tc_pres_wind.py').read())
exec(open(path+'interp_bdy.py').read())
exec(open(path+'solover.py').read())

# typhoon best track
f_tyBST = nc4.Dataset(tyBSTFile,'r')
Po     = f_tyBST.variables["Po"][:] 
Rm     = f_tyBST.variables["RMW"][:] 
latTy  = f_tyBST.variables["latTy"][:] 
lonTy  = f_tyBST.variables["lonTy"][:] 
VT     = f_tyBST.variables["VT"][:] 
Theta  = f_tyBST.variables["Theta"][:] 
timeTy = f_tyBST.variables["times"][:] 

for ID in range(dom_num):
    print("domain ",dom_id[ID])
    set_domain(domain_id=dom_id[ID])
    nt = int(run_hour*3600/dt[ID]) # max time steps
    ### origin location(0,0)
    Lon0_C = grid_C_Lon_Deg[0,0]; Lat0_C = grid_C_Lat_Deg[0,0]
    Lon0_U = grid_U_Lon_Deg[0,0]; Lat0_U = grid_U_Lat_Deg[0,0]
    Lon0_V = grid_V_Lon_Deg[0,0]; Lat0_V = grid_V_Lat_Deg[0,0]
    
    latTy0 = latTy[0]; lonTy0 = lonTy[0]
    Po0    = Po[0]
    Rm0    = Rm[0] #1000*1119.0 * (1010.0-950)**(-0.805) 
    VT0    = VT[0]
    Theta0 = Theta[0]
    VTx    = VT0*np.sin(Theta0*deg2rad)
    VTy    = VT0*np.cos(Theta0*deg2rad)
    
    print("typhoon location, lon, lat =",lonTy0,latTy0 )
    print("typhoon center pressure and Rmax, Po, Rm =",Po0,Rm0 )
    print("typhoon transfer velocity, VTx, VTy =",VTx,VTy )
    print("dt=%f" %(dt[ID]))
    print("dx=%f dy=%f" %(dx[ID],dy[ID]))
    print("dx max and min(km)=",np.max(dx[ID]*R*np.cos(grid_U_Lat)),np.min(dx[ID]*R*np.cos(grid_U_Lat)) )
    print("dy max and min(km)=",np.max(dy[ID]*R),np.min(dy[ID]*R) )
    
    presCC = get_ty_pres(gridXC,gridYC,Lon0_C,Lat0_C,lonTy0,latTy0,Po0,Rm0)
    windUU = get_ty_wind_stress(gridXU,gridYU,grid_U_Lat,Lon0_U,Lat0_U,VTx,VTy,lonTy0,latTy0,Po0,Rm0,output='windU')
    windVV = get_ty_wind_stress(gridXV,gridYV,grid_V_Lat,Lon0_V,Lat0_V,VTx,VTy,lonTy0,latTy0,Po0,Rm0,output='windV')
    
    # initial conditio
    u0   = np.zeros((Ny_u,Nx_u)) 
    v0   = np.zeros((Ny_v,Nx_v))
    #eta0 = np.zeros((Ny_c,Nx_c))
    eta0 = (Pe-presCC)/(g*rho_w)
    
    
    u0   = u0*maskU
    v0   = v0*maskV
    eta0 = eta0*maskC
    
    # save data
    ssm_out_file = "ssm_out_d0"+str(ID+1)+".nc"
    shutil.copyfile(gridFile[ID], ssm_out_file)
    fout = nc4.Dataset(ssm_out_file,mode='a')
    # creat time dimension
    fout.createDimension('times', None)
    u_save   = fout.createVariable('u','f8',('times','lat','lonStag'))
    v_save   = fout.createVariable('v','f8',('times','latStag','lon'))
    eta_save = fout.createVariable('eta','f8',('times','lat','lon'))
    fill_value = 9.99e+36 
    u_save.units   = "m*s-1"
    v_save.units   = "m*s-1"
    eta_save.units = "m"
    u_save.long_name   = "zonal velocity"
    v_save.long_name   = "meridional velocity"
    eta_save.long_name = "storm surge elevation"
    
    it_ty = 0
    ip  = 0
    latTy0 = latTy[it_ty]
    lonTy0 = lonTy[it_ty]
    latTy1 = latTy[it_ty+1]
    lonTy1 = lonTy[it_ty+1]
    Po0    = Po[it_ty]
    Po1    = Po[it_ty+1]
    Rm0    = Rm[it_ty]
    Rm1    = Rm[it_ty+1]
    VT0    = VT[it_ty]
    Theta0 = Theta[it_ty]
    VTx    = VT0*np.sin(Theta0*deg2rad)
    VTy    = VT0*np.cos(Theta0*deg2rad)
    lonTy01 = lonTy0 + (lonTy1-lonTy0)/360*ip
    latTy01 = latTy0 + (latTy1-latTy0)/360*ip
    Po01    = Po0 + (Po1-Po0)/360*ip
    Rm01    = Rm0 + (Rm1-Rm0)/360*ip
    
    spin_up = 1
    if spin_up == 1:
        print("spin up ...")
        nt_spin_up = int(6*3600/dt[ID])
        presC = get_ty_pres(gridXC,gridYC,Lon0_C,Lat0_C,lonTy0,latTy0,Po0,Rm0)
        windU = get_ty_wind_stress(gridXU,gridYU,grid_U_Lat,Lon0_U,Lat0_U,VTx,VTy,lonTy0,latTy0,Po0,Rm0,output='windU')
        windV = get_ty_wind_stress(gridXV,gridYV,grid_V_Lat,Lon0_V,Lat0_V,VTx,VTy,lonTy0,latTy0,Po0,Rm0,output='windV')
        if dom_id[ID] > 1:
            parent_file = "ssm_out_d01.nc"
            sub_file    = "ssm_out_d0"+str(dom_id[ID])+".nc"
            itBdy = 0   # initial time
            dictBdy = get_bdy_from_parent_domain(parent_file,sub_file,itBdy) 
        for it in range(nt_spin_up):
            ### u v eta forward
            u0,v0,eta0 = time_integration(u0,v0,eta0,presC,windU,windV,dx[ID],dy[ID],dt[ID],integration_method) 
            if dom_id[ID] == 1: 
                u0,v0,eta0 = LBC_Dom01(u0,v0,eta0,presC,uvLBC[ID],etaLBC[ID])
            else:
                u0,v0,eta0 = LBC_Dom02(u0,v0,eta0,dictBdy,uvLBC[ID],etaLBC[ID])
            itime = it*dt[ID]
        eta_max = np.max(np.abs(eta0))
        u_max   = np.max(np.abs(u0))
        v_max   = np.max(np.abs(v0))
        print("time = %6i seconds or %6.2f hours, ets = %4.6f u = %4.6f v = %4.6f lonTy = %3.3f latTy = %2.3f "  % (itime,itime/3600,eta_max,u_max,v_max,lonTy01,latTy01))
        print("end of spin up process")
    # save initial condition
    it_save = 0
    u_save0   = u0.copy()
    v_save0   = v0.copy()
    eta_save0 = eta0.copy()
    u_save0[maskU==0]   = fill_value
    v_save0[maskV==0]   = fill_value
    eta_save0[maskC==0] = fill_value
    u_save[it_save,:,:]   = u_save0
    v_save[it_save,:,:]   = v_save0
    eta_save[it_save,:,:] = eta_save0
     
    num_ty_interval = int(ty_interval*3600/dt[ID])
    itime = 0
    if dom_id[ID] > 1:
        parent_file = "ssm_out_d01.nc"
        sub_file    = "ssm_out_d0"+str(dom_id[ID])+".nc"
        itBdy = 0
        dictBdy = get_bdy_from_parent_domain(parent_file,sub_file,itBdy) 
    for it in range(1,nt+1):
        ip += 1
        if it%num_ty_interval == 0 :
            it_ty += 1
            it_ty2 = it_ty+1
            ip = 0 
            if it == nt:
                it_ty2 = it_ty
            latTy0 = latTy[it_ty]
            lonTy0 = lonTy[it_ty]
            latTy1 = latTy[it_ty2]
            lonTy1 = lonTy[it_ty2]
            Po0    = Po[it_ty]
            Po1    = Po[it_ty2]
            Rm0    = Rm[it_ty]
            Rm1    = Rm[it_ty2]
            VT0    = VT[it_ty]
            Theta0 = Theta[it_ty]
            VTx    = VT0*np.sin(Theta0*deg2rad)
            VTy    = VT0*np.cos(Theta0*deg2rad)
        lonTy01 = lonTy0 + (lonTy1-lonTy0)/num_ty_interval*ip
        latTy01 = latTy0 + (latTy1-latTy0)/num_ty_interval*ip
        Po01    = Po0 + (Po1-Po0)/num_ty_interval*ip
        Rm01    = Rm0 + (Rm1-Rm0)/num_ty_interval*ip
        presC = get_ty_pres(gridXC,gridYC,Lon0_C,Lat0_C,lonTy01,latTy01,Po01,Rm01)
        windU = get_ty_wind_stress(gridXU,gridYU,grid_U_Lat,Lon0_U,Lat0_U,VTx,VTy,lonTy01,latTy01,Po01,Rm01,output='windU')
        windV = get_ty_wind_stress(gridXV,gridYV,grid_V_Lat,Lon0_V,Lat0_V,VTx,VTy,lonTy01,latTy01,Po01,Rm01,output='windV')
     
        ### u v eta forward
        u0,v0,eta0 = time_integration(u0,v0,eta0,presC,windU,windV,dx[ID],dy[ID],dt[ID],integration_method)
        if dom_id[ID] == 1: 
            u0,v0,eta0 = LBC_Dom01(u0,v0,eta0,presC,uvLBC[ID],etaLBC[ID])
        else:
            if itime%(history_interval[ID]*60) == 0 and itime != 0:
                itBdy += 1  
                dictBdy = get_bdy_from_parent_domain(parent_file,sub_file,itBdy) 
            u0,v0,eta0 = LBC_Dom02(u0,v0,eta0,dictBdy,uvLBC[ID],etaLBC[ID])
        
        itime = it*dt[ID]
        # save data
        if itime%(history_interval[ID]*60) == 0:
            print("save data at",itime/3600,"hours")
            it_save +=1
            u_save0   = u0.copy()
            v_save0   = v0.copy()
            eta_save0 = eta0.copy()
            u_save0[maskU==0]   = fill_value
            v_save0[maskV==0]   = fill_value
            eta_save0[maskC==0] = fill_value
            u_save[it_save,:,:]   = u_save0
            v_save[it_save,:,:]   = v_save0
            eta_save[it_save,:,:] = eta_save0
            eta_max = np.max(np.abs(eta0))
            u_max   = np.max(np.abs(u0))
            v_max   = np.max(np.abs(v0))
            print("time = %6i seconds or %6.2f hours, ets = %4.6f u = %4.6f v = %4.6f lonTy = %3.3f latTy = %2.3f "  % (itime,itime/3600,eta_max,u_max,v_max,lonTy01,latTy01))
     
    fout.close()

# caculate comsume time
end_time = time.time()
run_time = end_time - start_time
print("run time = %6.2f seconds %6.3f minutes %6.4f hours" %(run_time,run_time/60,run_time/3600))


