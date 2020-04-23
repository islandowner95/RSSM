def time_integration(u_b,v_b,eta_b,pres_in,windU_in,windV_in,dx,dy,dt,method):
    if method == "RK3":
        u_a   = u_b.copy()
        v_a   = v_b.copy()
        eta_a = eta_b.copy()
        """
        # step 1
        u1,v1 = uv_forward(u_a,v_a,u_b,v_b,eta_b,pres_in,windU_in,windV_in,dt)
        eta1  = eta_forward(eta_a,u_b,v_b,eta_b,dt)
        # step 2
        u2,v2 = uv_forward(u_a,v_a,u1,v1,eta1,pres_in,windU_in,windV_in,dt)
        eta2  = eta_forward(eta_a,u1,v1,eta1,dt)
        u2   = 3.0/4.0*u_a   + 1.0/4.0*u2
        v2   = 3.0/4.0*v_a   + 1.0/4.0*v2
        eta2 = 3.0/4.0*eta_a + 1.0/4.0*eta2
        # step 3
        u_out,v_out = uv_forward(u_a,v_a,u2,v2,eta2,pres_in,windU_in,windV_in,dt=dt)
        eta_out     = eta_forward(eta_a,u2,v2,eta2,dt=dt)
        u_out   = 1.0/3.0*u_a   + 2.0/3.0*u_out
        v_out   = 1.0/3.0*v_a   + 2.0/3.0*v_out
        eta_out = 1.0/3.0*eta_a + 2.0/3.0*eta_out
        """
        u1,v1 = uv_forward(u_a,v_a,u_b,v_b,eta_b,pres_in,windU_in,windV_in,dx,dy,dt=dt/3.0)
        eta1  = eta_forward(eta_a,u_b,v_b,eta_b,dx,dy,dt=dt/3.0)
        u2,v2 = uv_forward(u_a,v_a,u1,v1,eta1,pres_in,windU_in,windV_in,dx,dy,dt=dt/2.0)
        eta2  = eta_forward(eta_a,u1,v1,eta1,dx,dy,dt=dt/2.0)
        u_out,v_out = uv_forward(u_a,v_a,u2,v2,eta2,pres_in,windU_in,windV_in,dx,dy,dt=dt)
        eta_out     = eta_forward(eta_a,u2,v2,eta2,dx,dy,dt=dt)

    elif method == "semi_implict":
        ### u v forward
        u_out,v_out = uv_forward(u_b,v_b,u_b,v_b,eta_b,pres_in,windU_in,windV_in,dx,dy,dt=dt)
        ### eta forward
        eta_out = eta_forward(eta_b,u_out,v_out,eta_b,dx,dy,dt=dt)
    else:
        print("Error: the method of integration is not exist")
        exit()
    return u_out,v_out,eta_out

def uv_forward(u_a,v_a,u_b,v_b,eta_b,pres_in,windU_in,windV_in,dx,dy,dt):
    """
    input:
         u_b[Ny_u,Nx_u]
         v_b[Ny_v,Nx_v]
         eta_b[Ny_c,Nx_c]
         pres_in[Ny_c,Nx_c]
    output:
         u_out[Ny_u,Nx_u]
         v_out[Ny_v,Nx_v]
         
    """ 
    u_out = u_a
    v_out = v_a
    h = eta_b + H # total water depth
    ### u tendency 
    dzdx  = g*(eta_b[1:-1,1:]-eta_b[1:-1,:-1])/(dx*R*np.cos(grid_U_Lat[1:-1,1:-1]))
    ududx = (u_b[1:-1,2:]**2 - u_b[1:-1,:-2]**2)/(dx*4.0*R*np.cos(grid_U_Lat[1:-1,1:-1]))  # d(1/2)u**2/dx
    vdudy = ((v_b[2:-1,:-1]+v_b[2:-1,1:])*(u_b[2:,1:-1]-u_b[1:-1,1:-1])+ \
             (v_b[1:-2,:-1]+v_b[1:-2,1:])*(u_b[1:-1,1:-1]-u_b[:-2,1:-1]))/(dy*R*4.0) 
    PGFu = (pres_in[1:-1,1:]-pres_in[1:-1,:-1])/(rho_w*dx*R*np.cos(grid_U_Lat[1:-1,1:-1]))
    #PGFu = 0
    ### v tendency
    dzdy  = g*(eta_b[1:,1:-1]-eta_b[:-1,1:-1])/(dy*R)
    udvdx = ((u_b[:-1,2:-1]+u_b[1:,2:-1])*(v_b[1:-1,2:]-v_b[1:-1,1:-1])+ \
             (u_b[:-1,1:-2]+u_b[1:,1:-2])*(v_b[1:-1,1:-1]-v_b[1:-1,:-2]))/ \
             (dx*4.0*R*np.cos(grid_V_Lat[1:-1,1:-1]))
    vdvdy = (v_b[2:,1:-1]**2 - v_b[:-2,1:-1]**2)/(dy*4.0*R) # d(1/2)v**2/dy
    PGFv = (pres_in[1:,1:-1]-pres_in[:-1,1:-1])/(rho_w*dy*R)
    #PGFv =0
    if coriolisForce == "explict":
        fv = f_U[1:-1,1:-1]*0.25*(v_b[1:-2,:-1]+v_b[1:-2,1:]+v_b[2:-1,:-1]+v_b[2:-1,1:]) 
        fu = f_V[1:-1,1:-1]*0.25*(u_b[:-1,1:-2]+u_b[:-1,2:-1]+u_b[1:,1:-2]+u_b[1:,2:-1])
        u_out[1:-1,1:-1] = u_out[1:-1,1:-1] + (-1.0*dzdx - ududx - vdudy + fv - PGFu) * dt
        v_out[1:-1,1:-1] = v_out[1:-1,1:-1] + (-1.0*dzdy - udvdx - vdvdy - fu - PGFv) * dt
    else:
        u_out[1:-1,1:-1] = u_out[1:-1,1:-1] + (-1.0*dzdx - ududx - vdudy - PGFu) * dt
        v_out[1:-1,1:-1] = v_out[1:-1,1:-1] + (-1.0*dzdy - udvdx - vdvdy - PGFv) * dt

    if windStress == 1:
        hu  = np.zeros((Ny_u,Nx_u)) 
        # u
        hu[:,1:-1] = (h[:,:-1]+h[:,1:])/2.0
        hu[:,[0,-1]] = hu[:,[1,-2]]
        hu[maskU==0] = 1
        u_out[:,:] += dt*(windU_in[:,:]/rho_w/hu)
        #taox = dt*windU_in[:,:]/rho_w/hu 
        # v
        hv  = np.zeros((Ny_v,Nx_v))
        hv[1:-1,:] = (h[:-1,:]+h[1:,:])/2.0
        hv[[0,-1],:] = hv[[1,-2],:]
        hv[maskV==0] = 1
        v_out[:,:] += dt*(windV_in[:,:]/rho_w/hv)

    if bottomDrag == 1:
        if bottomDrag_option ==1: # semi implict
            # u
            vU  = np.zeros((Ny_u,Nx_u)) # v interpolate to u position
            vU[:,1:-1] = (v_b[:-1,:-1]+v_b[:-1,1:]+v_b[1:,:-1]+v_b[1:,1:])/4.0 # v interpolate to u position
            vU[:,[0,-1]] = vU[:,[1,-2]]
            Rx = dragConst*dt*np.sqrt(u_b**2+vU**2)
            u_out = u_out/(1+Rx)
            # v
            uV  = np.zeros((Ny_v,Nx_v)) # u interpolate to v position
            uV[1:-1,:] = (u_b[:-1,:-1]+u_b[:-1,1:]+u_b[1:,:-1]+u_b[1:,1:])/4.0 # u interpolate to v position
            uV[[0,-1],:] = uV[[1,-2],:]
            Ry = dragConst*dt*np.sqrt(uV**2+v_b**2)
            v_out = v_out/(1+Ry)
        elif bottomDrag_option ==2: # explict
            # u
            vU  = np.zeros((Ny_u,Nx_u)) # v interpolate to u position
            vU[:,1:-1] = (v_b[:-1,:-1]+v_b[:-1,1:]+v_b[1:,:-1]+v_b[1:,1:])/4.0 # v interpolate to u position
            vU[:,[0,-1]] = vU[:,[1,-2]]
            uv_U = np.sqrt(u_b**2+vU**2) # total speed at U position
            u_frict = dragConst*u_b*uv_U
            # v
            uV  = np.zeros((Ny_v,Nx_v)) # u interpolate to v position
            uV[1:-1,:] = (u_b[:-1,:-1]+u_b[:-1,1:]+u_b[1:,:-1]+u_b[1:,1:])/4.0 # u interpolate to v position
            uV[[0,-1],:] = uV[[1,-2],:]
            uv_V = np.sqrt(v_b**2+uV**2) # total speed at V position
            v_frict = dragConst*v_b*uv_V
            u_out -= u_frict*dt
            v_out -= v_frict*dt
        else:
            print("Error: bottomDrag_option is not exits ")
            exit()
   
    if coriolisForce == "semi_implicit":
        vU  = np.zeros((Ny_u,Nx_u)) # v interpolate to u position
        vU[:,1:-1] = (v_b[:-1,:-1]+v_b[:-1,1:]+v_b[1:,:-1]+v_b[1:,1:])/4.0 # v interpolate to u position
        vU[:,[0,-1]] = vU[:,[1,-2]]
        uV  = np.zeros((Ny_v,Nx_v)) # u interpolate to v position
        uV[1:-1,:] = (u_b[:-1,:-1]+u_b[:-1,1:]+u_b[1:,:-1]+u_b[1:,1:])/4.0 # u interpolate to v position
        uV[[0,-1],:] = uV[[1,-2],:]
        alphaU = dt*f_U
        betaU  = 0.25*alphaU**2
        alphaV = dt*f_V
        betaV  = 0.25*alphaV**2
        u_out = (u_out-betaU*u_b+alphaU+alphaU*vU)/(1+betaU)
        v_out = (v_out-betaV*v_b+alphaV+alphaV*uV)/(1+betaV)
               
    if lateralDiffusion == 1 :
        heU = h[1:-1,1:]
        hwU = h[1:-1,:-1] 
        hnU = 0.25*(h[1:-1,:-1]+h[1:-1,1:]+h[2:,:-1]+h[2:,1:])
        hsU = 0.25*(h[:-2,:-1]+h[:-2,1:]+h[1:-1,:-1]+h[1:-1,1:])
        hU  = 0.50*(h[1:-1,:-1]+h[1:-1,1:])
        hU[maskU[1:-1,1:-1]==0] = 1 # Avoid division by zero
        diffU = Ah/hU*( ( heU*(u_b[1:-1,2:]-u_b[1:-1,1:-1])/(dx*R*np.cos(grid_U_Lat[1:-1,1:-1])) - \
                          hwU*(u_b[1:-1,1:-1]-u_b[1:-1,:-2])/(dx*R*np.cos(grid_U_Lat[1:-1,1:-1])) )/ \
                         (dx*R*np.cos(grid_U_Lat[1:-1,1:-1])) + \
                         ( hnU*(u_b[2:,1:-1]-u_b[1:-1,1:-1])/(dy*R) - \
                           hsU*(u_b[1:-1,1:-1]-u_b[:-2,1:-1])/(dy*R) )/ \
                         (dy*R) )
        u_out[1:-1,1:-1] += diffU*dt

        heV = 0.25*(h[:-1,1:-1]+h[:-1,2:]+h[1:,1:-1]+h[1:,2:])
        hwV = 0.25*(h[:-1,:-2]+h[:-1,1:-1]+h[1:,:-2]+h[1:,1:-1])
        hnV = h[1:,1:-1]
        hsV = h[:-1,1:-1]
        hV  = 0.50*(h[:-1,1:-1]+h[1:,1:-1])
        hV[maskV[1:-1,1:-1]==0] = 1 # Avoid division by zero
        diffV = Ah/hV*( ( heV*(v_b[1:-1,2:]-v_b[1:-1,1:-1])/(dx*R*np.cos(grid_V_Lat[1:-1,1:-1])) - \
                          hwV*(v_b[1:-1,1:-1]-v_b[1:-1,:-2])/(dx*R*np.cos(grid_V_Lat[1:-1,1:-1])) )/ \
                         (dx*R*np.cos(grid_V_Lat[1:-1,1:-1])) + \
                         ( hnV*(v_b[2:,1:-1]-v_b[1:-1,1:-1])/(dy*R) - \
                           hsV*(v_b[1:-1,1:-1]-v_b[:-2,1:-1])/ \
                         (dy*R) )/(dy*R) )
        v_out[1:-1,1:-1] += diffV*dt
    u_out = u_out*maskU
    v_out = v_out*maskV
    #print("dzdx=%f ududx=%1.12f vdudy=%1.12f PGFu=%f windU=%f" %( np.mean(np.abs(dzdx)), np.mean(np.abs(ududx)), np.mean(np.abs(vdudy)), np.mean(np.abs(PGFu)),np.mean(np.abs(taox)) ) )
    #print("max dzdy=%f udvdx=%1.12f vdvdy=%1.12f PGFv=%f windU=%f " %( np.max(np.abs(dzdy)), np.max(np.abs(udvdx)), np.max(np.abs(vdvdy)), np.max(np.abs(PGFv)), np.max(np.abs(taox)) ) )
    #print("dzdx=%f ududx=%1.12f vdudy=%1.12f PGFu=%f " %( np.mean(np.abs(dzdx)), np.mean(np.abs(ududx)), np.mean(np.abs(vdudy)), np.mean(np.abs(PGFu)) ) )
    #print("max dzdy=%f udvdx=%1.12f vdvdy=%1.12f PGFv=%f " %( np.max(np.abs(dzdy)), np.max(np.abs(udvdx)), np.max(np.abs(vdvdy)), np.max(np.abs(PGFv)) ) )
    return u_out, v_out


def eta_forward(eta_a,u_b,v_b,eta_b,dx,dy,dt):
    eta_out = eta_a
    h = eta_b + H # total water depth
    ### eta tendency
    """
    he  = np.zeros((Ny_c,Nx_c))
    hw  = np.zeros((Ny_c,Nx_c))
    hn  = np.zeros((Ny_c,Nx_c))
    hs  = np.zeros((Ny_c,Nx_c))
    he[:,:-1] = (h[:,1:]+h[:,:-1])/2.0
    he[:,-1]  = he[:,-2]
    hw[:,1:]  = (h[:,1:]+h[:,:-1])/2.0
    hw[:,0]   = hw[:,1]
    hn[:-1,:] = (h[1:,:]+h[:-1,:])/2.0
    hn[-1,:]  = hn[-2,:]
    hs[1:,:]  = (h[1:,:]+h[:-1,:])/2.0
    hs[0,:]   = hs[1,:]
    dudx = (u_b[:,1:]*he - u_b[:,:-1]*hw)/(dx*R*np.cos(grid_C_Lat))
    dvdy = (v_b[1:,:]*hn*np.cos(grid_V_Lat[1:,:]) - v_b[:-1,:]*hs*np.cos(grid_V_Lat[:-1,:]))/(dy*R*np.cos(grid_C_Lat))
    """
    he  = np.zeros((Ny_c,Nx_c))
    hw  = np.zeros((Ny_c,Nx_c))
    hn  = np.zeros((Ny_c,Nx_c))
    hs  = np.zeros((Ny_c,Nx_c))
    u_plus  = 0.5*(u_b+np.abs(u_b))
    u_minus = 0.5*(u_b-np.abs(u_b))
    v_plus  = 0.5*(v_b+np.abs(v_b))
    v_minus = 0.5*(v_b-np.abs(v_b))
    he[:,:-1] = h[:,1:];  he[:,-1] = h[:,-1] 
    hw[:,1:]  = h[:,:-1]; he[:,0]  = h[:,0] 
    hn[:-1,:] = h[1:,:];  hn[-1,:] = h[-1,:] 
    hs[1:,:]  = h[:-1,:]; hs[0,:]  = h[0,:] 
    dudx = (u_plus[:,1:]*h+u_minus[:,1:]*he - u_plus[:,:-1]*hw-u_minus[:,:-1]*h)/(dx*R*np.cos(grid_C_Lat))
    dvdy = (v_plus[1:,:]*h*np.cos(grid_V_Lat[1:,:])+v_minus[1:,:]*hn*np.cos(grid_V_Lat[1:,:]) - \
            v_plus[:-1,:]*hs*np.cos(grid_V_Lat[:-1,:])-v_minus[:-1,:]*h*np.cos(grid_V_Lat[:-1,:]))/(dy*R*np.cos(grid_C_Lat))
    #print(np.max(np.abs(dudx)),np.max(np.abs(dvdy)))
    eta_out = eta_out - dt*(dudx+dvdy)
    eta_out = eta_out*maskC
    #print("dudx=%1.12f dvdy=%1.12f" %( np.mean(np.abs(dudx)), np.mean(np.abs(dvdy)) ) )
    #print("max dudx=%1.12f dvdy=%1.12f" %( np.max(np.abs(dudx)), np.max(np.abs(dvdy)) ) )
    # 2D 1-order shapiro filte
    if shapiroFilte == 1:
        eps = 0.5
        eta_out[1:-1,1:-1] = (1-eps)*eta_out[1:-1,1:-1]+ \
                      0.25*eps*(eta_out[0:-2,1:-1]+eta_out[2:,1:-1]+ \
                                eta_out[1:-1,0:-2]+eta_out[1:-1,2:])
    return eta_out




def LBC_Dom01(u,v,eta,pres,uvLBC,etaLBC):        
    # u v boundary condition
    if uvLBC == 1:
        u[:,[0,-1]] = 0.0 # u[:,[-1,2]]
        v[:,[0,-1]] = 0.0 #v[:,[-1,2]]
        u[[0,-1],:] = 0.0
        v[[0,-1],:] = 0.0
    elif uvLBC == 2:
        u[:,[0,-1]] = 2.0*u[:,[1,-2]]-u[:,[2,-3]]
        v[:,[0,-1]] = 2.0*v[:,[1,-2]]-v[:,[2,-3]]
        u[[0,-1],:] = 2.0*u[[1,-2],:]-u[[2,-3],:]
        v[[0,-1],:] = 2.0*v[[1,-2],:]-v[[2,-3],:]
    else:
        print("Error: uvLBC option ")
        exit()
    
    if etaLBC == 1:
        eta[:,[0,-1]] = 0.0 
        eta[[0,-1],:] = 0.0
    elif etaLBC == 2:
        eta_pres = (Pe-pres)/(g*rho_w)
        eta[[0,1,2,-1,-2,-3],:] = eta_pres[[0,1,2,-1,-2,-3],:]
        eta[:,[0,1,2,-1,-2,-3]] = eta_pres[:,[0,1,2,-1,-2,-3]]
    else:
        print("Error: etaLBC option ")
        exit()
    u   = u*maskU
    v   = v*maskV
    eta = eta*maskC
    return u, v, eta


def LBC_Dom02(u,v,eta,dictBdy,uvLBC,etaLBC):        
    # u v boundary condition
    eta_N = dictBdy["eta_N"] 
    eta_S = dictBdy["eta_S"] 
    eta_E = dictBdy["eta_E"] 
    eta_W = dictBdy["eta_W"] 
    u_N = dictBdy["u_N"] 
    u_S = dictBdy["u_S"] 
    u_E = dictBdy["u_E"] 
    u_W = dictBdy["u_W"] 
    v_N = dictBdy["v_N"] 
    v_S = dictBdy["v_S"] 
    v_E = dictBdy["v_E"] 
    v_W = dictBdy["v_W"] 
    if uvLBC == 1:
        u[-1,:] = u_N
        u[-2,:] = 0.75*u_N + 0.25*u[-4,:]
        u[-3,:] = 0.25*u_N + 0.75*u[-4,:]
        u[0,:]  = u_S
        u[1,:]  = 0.75*u_S + 0.25*u[3,:]
        u[2,:]  = 0.25*u_S + 0.75*u[3,:]
        u[:,-1] = u_E
        u[:,-2] = 0.75*u_E + 0.25*u[:,-4]
        u[:,-3] = 0.25*u_E + 0.75*u[:,-4]
        u[:,0]  = u_W
        u[:,1]  = 0.75*u_W + 0.25*u[:,3]
        u[:,2]  = 0.25*u_W + 0.75*u[:,3]

        v[-1,:] = v_N
        v[-2,:] = 0.75*v_N + 0.25*v[-4,:]
        v[-3,:] = 0.25*v_N + 0.75*v[-4,:]
        v[0,:]  = v_S
        v[1,:]  = 0.75*v_S + 0.25*v[3,:]
        v[2,:]  = 0.25*v_S + 0.75*v[3,:]
        v[:,-1] = v_E
        v[:,-2] = 0.75*v_E + 0.25*v[:,-4]
        v[:,-3] = 0.25*v_E + 0.75*v[:,-4]
        v[:,0]  = v_W
        v[:,1]  = 0.75*v_W + 0.25*v[:,3]
        v[:,2]  = 0.25*v_W + 0.75*v[:,3]
    else:
        print("Error: uvLBC option ")
        exit()
    if etaLBC == 1:
        eta[-1,:] = eta_N
        eta[-2,:] = 0.75*eta_N + 0.25*eta[-4,:]
        eta[-3,:] = 0.25*eta_N + 0.75*eta[-4,:]
        eta[0,:]  = eta_S
        eta[1,:]  = 0.75*eta_S + 0.25*eta[3,:]
        eta[2,:]  = 0.25*eta_S + 0.75*eta[3,:]
        eta[:,-1] = eta_E
        eta[:,-2] = 0.75*eta_E + 0.25*eta[:,-4]
        eta[:,-3] = 0.25*eta_E + 0.75*eta[:,-4]
        eta[:,0]  = eta_W
        eta[:,1]  = 0.75*eta_W + 0.25*eta[:,3]
        eta[:,2]  = 0.25*eta_W + 0.75*eta[:,3]
    else:
        print("Error: etaLBC option ")
        exit()
    u   = u*maskU
    v   = v*maskV
    eta = eta*maskC
    return u, v, eta



