import numpy as np
from namelist_python import read_namelist_file
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpathes
import matplotlib.patches as patches
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon


# This function is to plot the base map
def plotBase(fig, latMin,latMax,lonMin,lonMax):
    m = Basemap(projection='merc',
                lon_0=0,lat_0=0,lat_ts=0,
                llcrnrlat=latMin,urcrnrlat=latMax,
                llcrnrlon=lonMin,urcrnrlon=lonMax,
                resolution='c')
    #m.drawcountries(linewidth=1, color='k')
    #m.drawmapscale(-90, 5, -90, 5, 1000, barstyle='fancy')
    m.bluemarble(scale=1)
    return m

def set_lonlat(_m, lon_list, lat_list, lon_labels, lat_labels, lonlat_size):
    """
    :param _m: Basemap实例
    :param lon_list: 经度 详见Basemap.drawmeridians函数>介绍
    :param lat_list: 纬度 同上
    :param lon_labels: 标注位置 [左, 右, 上, 下] bool值 默认只标注左上待完善 可使用twinx和twiny实现
    :param lat_labels: 同上
    :param lonlat_size: 字体大小
    :return:
    """
    lon_dict = _m.drawmeridians(lon_list, labels=lon_labels, color='none', fontsize=lonlat_size)
    lat_dict = _m.drawparallels(lat_list, labels=lat_labels, color='none', fontsize=lonlat_size)
    lon_list = []
    lat_list = []
    for lon_key in lon_dict.keys():
        try:
            lon_list.append(lon_dict[lon_key][1][0].get_position()[0])
        except:
            continue

    for lat_key in lat_dict.keys():
        try:
            lat_list.append(lat_dict[lat_key][1][0].get_position()[1])
        except:
            continue
    ax = plt.gca()
    ax.xaxis.tick_top()
    ax.set_yticks(lat_list)
    ax.set_xticks(lon_list)
    ax.tick_params(labelcolor='none')
    return None

def draw_domain_rect(m, latMin, lonMin, latMax, lonMax, text):
    x0, y0 = m(lonMin, latMin)
    x1, y1 = m(lonMax,latMax)
    width  = x1-x0
    height = y1-y0
    rect = patches.Rectangle((x0,y0),width,height,linewidth=1,edgecolor='r',facecolor='none')
    ax = plt.gca()
    ax.add_patch(rect)
    xt, yt = m(lonMin+0.8, latMax-0.8)
    plt.text(xt,yt,text,color="white")
    return None

path = '../'
exec(open(path+'constants.py').read())
get_constants()
exec(open(path+'read_namelist.py').read())
read_namelist()


fig = plt.figure()
ax = fig.add_subplot()

m = plotBase(ax, latMin[0]-1,latMax[0]+1,lonMin[0]-1,lonMax[0]+1)
parallels = np.arange(latMin[0],latMax[0]+1,2)
meridians = np.arange(lonMin[0],lonMax[0]+1,2)
set_lonlat(m,meridians,parallels,[0,0,0,1], [1,0,0,0],10)

for i in range(len(latMax)):
    text = "domain "+str(i+1)
    draw_domain_rect(m, latMin[i], lonMin[i],latMax[i],lonMax[i], text)
ax.set_title("simulation domain",fontsize=12)
plt.show()
figName = "domain.png"
fig.savefig(figName)
plt.close()








