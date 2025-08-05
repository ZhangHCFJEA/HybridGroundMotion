# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 08:29:07 2021

@author: ZhangHC
"""

import os
from obspy.core import read
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

xytmpfile = 'surface_st.txt'
path = 'E:/kiknet_ml/surface'
files = os.listdir(path)
info, Lat, Lon = [], [], []
for f in range(0, len(files), 3):
    print(f)
    tr = read(os.path.join(path, files[f]))
    epcLat, epcLon, mag, dep = tr[0].stats.knet.evla, tr[0].stats.knet.evlo, tr[0].stats.knet.mag, tr[0].stats.knet.evdp
    stLat, stLon = tr[0].stats.knet.stla, tr[0].stats.knet.stlo
    #
    Lat.append(stLat)
    Lon.append(stLon)
    #
    # tr[0].detrend()
    # tr[0].stats.sampling_rate
    # accmax = tr[0].stats.knet.accmax #tr[0].data.max()*tr[0].stats.calib
    # info.append([epcLat, epcLon, mag, dep, stLat, stLon, accmax])
    

np.savetxt(os.path.join( 'E:/kiknet_ml/', xytmpfile), np.c_[Lon, Lat])

#    
np.savetxt('E:/kiknet_ml/kik-info.txt', info)

minLon, minLat, maxLon, maxLat = 128, 29, 150, 47
fig = plt.figure(1, figsize=(10, 10))
map = Basemap(llcrnrlon=minLon, llcrnrlat=minLat, urcrnrlon=maxLon, urcrnrlat=maxLat, 
              projection='merc', lat_0=39, lon_0=138, resolution='f')
map.drawparallels(np.arange(minLat, maxLat, 5),labels=[1,0,0,1])
map.drawmeridians(np.arange(minLon, maxLon, 5),labels=[1,0,0,1])
# map.etopo()
map.shadedrelief()
map.drawcoastlines()
for e in range(len(info)):
    x, y = map(info[e][1], info[e][0])
    dep, mag = info[e][3], info[e][2]
    # map.scatter(x, y, marker='o', s=mag*0.1, c=dep, cmap='hsv') 
    map.scatter(x, y, marker='o', color=dep, cmap='hsv') 

map.colorbar(label='focal depth (km)')
map.drawmapscale(147, 31, 0.1, 39.5, 500, barstyle='fancy')
plt.savefig('E:/kiknet_ml/kik-net.png', dpi=200, bbox_inches='tight')
plt.savefig('E:/kiknet_ml/kik-net.eps', dpi=200, bbox_inches='tight')
plt.show()


# map.colorbar(label='focal depth (m)')
# map.drawmapscale(-119.2, 32.5, 0.25, 39.5, 100, barstyle='fancy')
# map.drawmapscale(-117.8, 32.5, 4.25, 39.5, 100, fontsize = 10)
# x, y = map(-118.15, 34.14)
# plt.text(x, y, 'LA', fontsize=16, fontweight='bold', ha='center',va='center')

    
plt.plot(info[:][1], info[:][0], 'r.', markersize=2)

for e in range(len(info)):
    plt.plot(info[e][1], info[e][0], 'r.')
    