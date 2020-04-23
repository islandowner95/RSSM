import numpy as np
import math
from math import sin,radians,cos,asin,sqrt,atan2
import datetime
import time
from namelist_python import read_namelist_file

def getDistFromTy(lonTy,latTy,gridX,gridY):
    """
    get distance between typhoon center and all points
    output: distance (km)
    """
    EARTH_REDIUS = 6378.135
    latTyRad = deg2rad(latTy)
    gridYRad = deg2rad(gridY)
    a = latTyRad - gridYRad
    b = deg2rad(lonTy) - deg2rad(gridX)
    arg1 = np.sqrt(pow(np.sin(a/2), 2) + np.cos(latTyRad) * np.cos(gridYRad) * pow(np.sin(b/2), 2))
    m,n = np.shape(arg1)
    arg2 = np.zeros([m,n])
    for i in range(m):
        for j in range(n):
            arg2[i,j] = math.asin(arg1[i,j]) 
    dist =2*arg2 * EARTH_REDIUS
    return dist

def SphereDistance(lon1, lat1, lon2, lat2, option):
    """
    option: "deg" or "rad"
    """
    # This function compute the great-circle distance 
    # between (lat1, lon1) and (lat2, lon2) on
    # a sphere with given radius.
    radius = 6378.135E3 # radius of Earth, unit:m
    # degree to radians
    if option =="deg":
        lon1, lat1,lon2, lat2 = map(radians,[lon1, lat1,lon2, lat2])
    elif option =="rad":
        pass
    else:
        print("Error: option must be deg or rad.")
        exit()
    dlon = lon2 -lon1
    dlat = lat2 -lat1
    arg  = sin(dlat*0.5)**2 + cos(lat1)*cos(lat2)*sin(dlon*0.5)**2
    dist = 2.0 * radius * asin(sqrt(arg))
    return dist

def getAzimuth(lon1, lat1, lon2, lat2):
    # This function compute the azimuth from A(lat1, lon1) to
    # B(lat2, lon2) on lats and lons on degree
    lon1, lat1, lon2, lat2 = map(radians,[lon1,lat1,lon2,lat2])
    cosC = cos(90-lat2)*cos(90-lat1) +  \
         sin(90-lat2)*sin(90-lat1)*cos(lon2-lon1)
    sinC = sqrt(1-cosC*cosC)
    if sinC == 0:
          print(lon1, lat1, lon2, lat2)
          sinC = 0.0001
    arg1 = (sin(90-lat2)*sin(lon2-lon1))/sinC
    A = asin(arg1)
    A = A*180/math.pi
    if lat2 >= lat1:
        if lon2 >= lon1:
            Azimuth = A
        else:
            Azimuth = 360 + A
    else:
        Azimuth = 180 - A
    return Azimuth

def deg2rad(d):
    """
    degree to radian
    """
    return d * math.pi / 180.0

def get_time_dimension(startDate,endDate,interval):
    """
    input: starDate: '2018-07-08_00:00:00'
           interval: unit:hour
    """
    startTime = datetime.datetime.strptime(startDate, "%Y-%m-%d_%H:%M:%S")
    endTime   = datetime.datetime.strptime(endDate, "%Y-%m-%d_%H:%M:%S")
    total_seconds = (endTime - startTime).total_seconds()
    time = int(total_seconds/3600.0/interval)+1 
    return time

def get_total_time(startDate,endDate):
    startTime = datetime.datetime.strptime(startDate, "%Y-%m-%d_%H:%M:%S")
    endTime   = datetime.datetime.strptime(endDate, "%Y-%m-%d_%H:%M:%S")
    total_seconds = (endTime - startTime).total_seconds()
    return total_seconds

def read_typhoon_record():
    """
    read CMA-SIT typhoon bset track data
    output:
        presTy: pressure at typhoon center, unit:Pa
        RMW: typhoon maximum wind radius, unit:m
        lat and lon
    """
    namelist = read_namelist_file('namelist.ssm')
    file_bst = namelist.groups["typhoon"]["BSTdata"]
    ty_interval = namelist.groups["typhoon"]["ty_interval"]
    ty_interval_s = ty_interval*3600 # hour to second
    dateTy = [] # 2011060600
    rmwTy  = [] # km
    latTy  = []
    lonTy  = [] 
    presTy = [] # hPa
    VT     = [] # m/s
    Theta  = [] # degree
    f = open(file_bst,'r')
    lines = f.readlines()
    f.close()
    for line in lines:
        lineSplit = line.split()
        # 2018060700 -> 2018-06-07_00:00:00
        dateRec = lineSplit[0]
        dateRec = time.strptime(dateRec, "%Y%m%d%H")
        dateRec = time.strftime("%Y-%m-%d_%H:%M:%S", dateRec)
        latRec  = float(lineSplit[2])*0.1
        lonRec  = float(lineSplit[3])*0.1
        presRec = float(lineSplit[4])
        RmaxRec1 = 1119.0 * (1010.0-presRec)**(-0.805)
        RmaxRec2 = 28.52*np.tanh(0.0873*(latRec-28)) + 12.22*np.exp((presRec-1013.2)/33.86) + 0.2*15+37.22
        RmaxRec3 = np.exp(2.636-(5.086E-5)*(1010.0-presRec)**2+0.0394899*latRec) 
        RmaxRec4 = np.exp( 2.0633+0.0182*(1010.0-presRec)-1.9008E-4*(1010.0-presRec)**2+7.336E-4*latRec**2 ) 
        print(RmaxRec1,RmaxRec2,RmaxRec3,RmaxRec4)
        RmaxRec =  RmaxRec1
        dateTy.append(dateRec)
        latTy.append(latRec)
        lonTy.append(lonRec)
        presTy.append(presRec)
        rmwTy.append(RmaxRec)
    latTy  = np.array(latTy)
    lonTy  = np.array(lonTy)
    presTy = np.array(presTy)
    rmwTy  = np.array(rmwTy)
    presTy = presTy*100 # hPa -> Pa
    rmwTy  = rmwTy*1000  # km -> m
    n = np.shape(latTy)[0]
    VT    = np.zeros(n)
    Theta = np.zeros(n)
    for i in range(0,n-1):
        dist     = SphereDistance(lonTy[i],latTy[i],lonTy[i+1],latTy[i+1],option="deg")
        VT[i]    = dist/ty_interval_s
        Theta[i] = getAzimuth(lonTy[i],latTy[i],lonTy[i+1],latTy[i+1])
    VT[-1] = VT[-2]
    Theta[-1] = Theta[-2]
    dictTy = {}
    dictTy["dateTy"] = dateTy
    dictTy["rmwTy"]  = rmwTy
    dictTy["latTy"]  = latTy
    dictTy["lonTy"]  = lonTy
    dictTy["presTy"] = presTy
    dictTy["VT"]     = VT
    dictTy["Theta"]  = Theta

    return dictTy


if __name__ == '__main__':
    print("tools")
    #lonTy = 123.0
    #latTy = 30.0  
    lonTy = 133.6
    latTy = 21.8



    lonMin1 = 119.0
    lonMax1 = 127.0
    latMin1 = 25.0
    latMax1 = 35.0
    lonMin2 = 120.0
    lonMax2 = 125.0
    latMin2 = 27.0
    latMax2 = 33.0
    dx1     = 1/10.0
    dx2     = 1/30.0


    #print(np.min(distXU1))
    #print(np.max(distXU1))
    from namelist_python import read_namelist_file
    namelist = read_namelist_file('namelist.ssm')
    startDate = namelist.groups["share"]["startDate"]
    endDate   = namelist.groups["share"]["endDate"]
    interval  = namelist.groups["share"]["interval"]
    time =  get_time_dimension(startDate,endDate,interval)
    print(time)

    Theta = getAzimuth(120, 30, 120, 31)
    print(Theta)
    Theta = getAzimuth(120, 30, 121, 30)
    print(Theta)

