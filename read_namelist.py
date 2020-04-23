def read_namelist():
    global share,domain,typhoon,dynamics
    namelist  = read_namelist_file('namelist.ssm')
    share     = namelist.groups["share"]
    domain    = namelist.groups["domain"]
    typhoon   = namelist.groups["typhoon"]
    dynamics  = namelist.groups["dynamics"]
    read_share()
    read_domain()
    read_typhoon()
    read_dynamics()
    return None

def read_share():
    global startDate,endDate,interval,integration_method,run_hour
    global bottomDrag_option
    startDate = share["startDate"]
    endDate   = share["endDate"]
    interval  = share["interval"]
    run_hour  = share["run_hour"]
    integration_method = share["integration_method"]
    if integration_method == "semi_implict":
        bottomDrag_option = 1
    elif integration_method == "RK3":
        bottomDrag_option = 2
    else:
        print("Error: integration methid is not exsit")
        exit()
    return None

def read_domain():
    global dom_num,dom_id
    global latMin,latMax,lonMin,lonMax,dx,dy,dt
    global etopoFile,gridFile
    global uvLBC,etaLBC
    global history_interval
    dom_num = domain["dom_num"]
    if dom_num >1:
        dom_id = domain["dom_id"]
        latMin = domain["latMin"]
        latMax = domain["latMax"]
        lonMin = domain["lonMin"]
        lonMax = domain["lonMax"]
        dx = domain["dx"]
        dy = domain["dy"]
        dx = list(np.array(dx)*deg2rad)
        dy = list(np.array(dy)*deg2rad)
        dt = domain["dt"]
        gridFile = domain["gridFile"]
        uvLBC  = domain["uvLBC"]
        etaLBC = domain["etaLBC"]
        history_interval = domain["history_interval"]
    else: # transfer to list
        dom_id = []
        latMin = []; latMax = []
        lonMin = []; lonMax = []
        dx = []; dy = []; dt = []
        gridFile = []
        uvLBC = []; etaLBC = []
        history_interval = []
        dom_id.append(domain["dom_id"])
        latMin.append(domain["latMin"])
        latMax.append(domain["latMax"])
        lonMin.append(domain["lonMin"])
        lonMax.append(domain["lonMax"])
        dx.append(domain["dx"]*deg2rad)
        dy.append(domain["dy"]*deg2rad)
        dt.append(domain["dt"])
        gridFile.append(domain["gridFile"])
        uvLBC.append(domain["uvLBC"])
        etaLBC.append(domain["etaLBC"])
        history_interval.append(domain["history_interval"]) # minutes
    etopoFile = domain["etopoFile"]
    
    return None

def read_typhoon():
    global BSTdata,ty_interval,tyBSTFile,ty_pres_wind
    BSTdata      = typhoon["BSTdata"]
    ty_interval  = typhoon["ty_interval"]
    ty_pres_wind = typhoon["ty_pres_wind"]
    tyBSTFile    = typhoon["tyBSTFile"]
    return None

def read_dynamics():
    global windStress,bottomDrag,lateralDiffusion,shapiroFilte,coriolisForce
    windStress       = dynamics["windStress"]
    bottomDrag       = dynamics["bottomDrag"]
    lateralDiffusion = dynamics["lateralDiffusion"]
    shapiroFilte     = dynamics["shapiroFilte"]
    coriolisForce    = dynamics["coriolisForce"]
    if coriolisForce == "explict":
        pass
    elif coriolisForce == "semi_implict":
        pass
    else:
        print("Error: corilisForce = "+coriolisForce+" is not exsit")
        exit()
    return None


