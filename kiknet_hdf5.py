# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 16:42:25 2021

@author: ZhangHC
"""

import os, sys
sys.path.append(r'E:/python')
import JMATravelTime
from obspy.core import read, UTCDateTime
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from obspy.signal.trigger import ar_pick
from obspy.taup import TauPyModel
from obspy.geodetics.base import gps2dist_azimuth as g2d
from obspy.geodetics.base import kilometer2degrees as k2d
from scipy.signal import butter, filtfilt
import h5py
import subprocess
from shlex import split
import shutil
from gmpe_tools import bssa14_one_station, JMA_II, MMI_Worden2012, WGRW12, save_h5, tau_c
import pandas as pd

#
vs30_grdfile = 'E:/GMTtest/global_vs30_grd/global_vs30.grd'
xy_vs30tmpfile = 'vs30.xyz'
# xytmpfile = 'tmp.txt'
# np.savetxt(xytmpfile, np.c_[stLon, stLat])
stList = 'G:/kiknet_ml/surface_st.txt'
command = split('grdtrack ' + stList + ' -G' + vs30_grdfile + ' > ' + xy_vs30tmpfile)
p = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
res,err = p.communicate()
Vs = np.genfromtxt(xy_vs30tmpfile, usecols=[2])

# 
# path = 'G:/kiknet_ml/borehole'
path = 'F:/kiknet_ml/surface'
savePath = 'F:/kiknet_ml'
hdfFileName = 'kiknet_new.hdf5'
files = os.listdir(path)
model = TauPyModel(model="ak135")
f1, f2 = 1, 20
lta_p, sta_p, lta_s, sta_s = 1.0, 0.1, 4.0, 1.0
m_p, m_s, l_p, l_s = 2, 8, 0.1, 0.2
c = 1
hpf = 0.075
smp = 100
L = 10

# 
for f in range(60000, len(files), 3):
    print(f)
    tr = read(os.path.join(path, files[f]))
    tr += read(os.path.join(path, files[f+1]))
    tr += read(os.path.join(path, files[f+2]))
    df = tr[0].stats.sampling_rate
    #
    sta, net = tr[0].stats.station, tr[0].stats.network
    epcLat, epcLon, mag, dep = tr[0].stats.knet.evla, tr[0].stats.knet.evlo, tr[0].stats.knet.mag, tr[0].stats.knet.evdp
    stLat, stLon, stElv = tr[0].stats.knet.stla, tr[0].stats.knet.stlo, tr[0].stats.knet.stel
    startTime = tr[0].stats.starttime
    #
    index = f//3
    vs30 = Vs[index]

    #
    dist = g2d(epcLat, epcLon, stLat, stLon)
    hypoDist = ((dist[0]/1000.)**2+(dep**2))**0.5
    pga_bssa14 = bssa14_one_station(mag, hypoDist, vs30, intensity_measure='PGA') # g
    pgv_bssa14 = bssa14_one_station(mag, hypoDist, vs30, intensity_measure='PGV') # m/s
    
    
    # arrivals = model.get_travel_times(source_depth_in_km = dep,
    #                               distance_in_degree = k2d(dist[0]/1000.) )
    #
    tr.detrend()

    #
    tr[0].data = tr[0].data*tr[0].stats.calib
    tr[1].data = tr[1].data*tr[1].stats.calib
    tr[2].data = tr[2].data*tr[2].stats.calib
    
    year = tr[0].stats.starttime.year
    lat0 = tr[0].stats.knet.stla
    lon0 = tr[0].stats.knet.stlo
    lat1 = tr[0].stats.knet.evla
    lon1 = tr[0].stats.knet.evlo
    
    if year < 2020 :
        if startTime < UTCDateTime('2011-03-11T00:00') or startTime > UTCDateTime('2011-03-21T00:00'):
            try :
                filePath = 'D:/桌面/数据/JMA/h2000/h'+str(year)+'.csv'
                cataList = pd.read_csv(filePath, header=None)
                    
                tmp = []
                for s in range(len(cataList)):
                    tmp.append(np.abs(UTCDateTime(cataList[10][s])-tr[0].stats.knet.evot-9*60*60))
                    
                # multiple?
                k1 = np.where(tmp == np.min(tmp))[0]
                for k in range(len(k1)):
                    num = k1[k]
                    lat2 = cataList[6][num]
                    lon2 = cataList[7][num]
                    dist = g2d(lat1, lon1, lat2, lon2)
                        
                    if dist[0]/1000. < 10:
                        originTime = UTCDateTime(cataList[10][num])-9*60*60
                        model = 'tjma2001'
                        stDist = g2d(lat0, lon0, lat2, lon2)
                        travelTime = JMATravelTime.JMATravelTime(model, tr[0].stats.knet.evdp, stDist[0]/1000.)
                        timeP = travelTime[0]
                        break            
                
                # p_pick, s_pick = ar_pick(tr[0].data, tr[1].data, tr[2].data, df, f1, f2, lta_p, sta_p, 
                #                          lta_s, sta_s, m_p, m_s, l_p, l_s) #s_pick=True
                # p_theory = tr[0].stats.knet.evot+arrivals[0].time-tr[0].stats.starttime
                p_pick = originTime+timeP-tr[0].stats.starttime
                
                #
                pp = int(p_pick*df)
                if pp > 0 :
                    snr = np.std(tr[2].data[pp:int(pp+1*df)])/np.std(tr[2].data[int(pp-1*df):pp])
                    
                    JMA_inten = JMA_II(tr[0].data, tr[1].data, tr[2].data, df)
                    # plt.plot(tr[0].data)
                    # plt.plot([pp, pp], [-0.15, 0.15], color='r', linestyle='--')
                    # plt.xlim(pp-3*tr[0].stats.sampling_rate, pp+3*3*tr[0].stats.sampling_rate)
                    
                    
                    if df > smp:
                        tr.resample(smp)
            
                    # if snr > 1:
                    #
                    trv = tr.copy()
                    trv = trv.integrate()
                    trv.filter('highpass', freq=hpf, corners=4, zerophase=True)
                    
                    trd = trv.copy()
                    trd = trd.integrate()
                    trd.filter('highpass', freq=hpf, corners=4, zerophase=True)
                    
                    #
                    pga = np.max([tr[0].stats.knet.accmax, tr[1].stats.knet.accmax, tr[2].stats.knet.accmax])/1e2 #m/s/s
                    pgv = np.max([trv[0].data.max(), trv[1].data.max(), trv[2].data.max(),
                                  abs(trv[0].data.min()), abs(trv[0].data.min()), abs(trv[0].data.min())]) #m/s
                    # MMI
                    mmi_pgv = MMI_Worden2012(pgv*100., hypoDist, float(mag), 'PGV')
                    mmi_pga = MMI_Worden2012(pga*100., hypoDist, float(mag), 'PGA')
                    mmi_pgv2 = WGRW12(pgv*100., 1)
                    mmi_pga2 = WGRW12(pga*100., 0)
                    
                    # Tau_c & Pd
                    ppp = int(p_pick*smp) # resampled
                    disp = trd[2].data
                    vel = trv[2].data
                    tw = 3
                    Pd = max(abs(disp[ppp:int(ppp+tw*smp)]))
                    tc = tau_c(disp[ppp:int(ppp+tw*smp)], vel[ppp:int(ppp+tw*smp)])
                        
                    #
                    A = np.vstack((tr[0].data[ppp:int(ppp+L*smp)], tr[1].data[ppp:int(ppp+L*smp)], tr[2].data[ppp:int(ppp+L*smp)]))
                    V = np.vstack((trv[0].data[ppp:int(ppp+L*smp)], trv[1].data[ppp:int(ppp+L*smp)], trv[2].data[ppp:int(ppp+L*smp)]))
                    #
                    hf = h5py.File(os.path.join(savePath, hdfFileName), 'a')
                        
                    save_h5(hf, data=np.array([A]), target='accData')
                    save_h5(hf, data=np.array([V]), target='velData')
                    # save_h5(hf, data=np.array([eventid.encode('utf-8')]), target='eventID')
                    save_h5(hf, data=np.array([files[f].encode('utf-8')]), target='recordID')
                    save_h5(hf, data=np.array([float(mag)]), target='Mag')
                    save_h5(hf, data=np.array([float(pga)]), target='PGA') 
                    save_h5(hf, data=np.array([float(pgv)]), target='PGV') 
                    save_h5(hf, data=np.array([float(pga_bssa14[0]*10)]), target='PGA_BSSA14') 
                    save_h5(hf, data=np.array([float(pgv_bssa14[0])]), target='PGV_BSSA14')
                    save_h5(hf, data=np.array([float(mmi_pga)]), target='MMI_PGA') 
                    save_h5(hf, data=np.array([float(mmi_pgv)]), target='MMI_PGV') 
                    save_h5(hf, data=np.array([float(mmi_pga2)]), target='MMI_PGA2') 
                    save_h5(hf, data=np.array([float(mmi_pgv2)]), target='MMI_PGV2') 
                    save_h5(hf, data=np.array([float(JMA_inten)]), target='JMA intensity') 
                    save_h5(hf, data=np.array([float(vs30)]), target='Vs30') 
                    save_h5(hf, data=np.array([float(Pd)]), target='Pd') 
                    save_h5(hf, data=np.array([float(tc)]), target='Tau_c')
                    # save_h5(hf, data=np.array([tp_max]), target='Tau_p_max')
                    save_h5(hf, data=np.array([float(dist[0]/1000.)]), target='epicentral distance') 
                    save_h5(hf, data=np.array([float(dep)]), target='focal depth') 
                    # save_h5(hf1, data=np.array([arrival-st2[0].stats.starttime-60]), target='arrival') 
                    save_h5(hf, data=np.array([int(smp)]), target='sampling rate') 
                    save_h5(hf, data=np.array([sta.encode('utf-8')]), target='station') 
                    save_h5(hf, data=np.array([net.encode('utf-8')]), target='network')
                    #
                    save_h5(hf, data=np.array([float(epcLat)]), target='event_latitude')
                    save_h5(hf, data=np.array([float(epcLon)]), target='event_longitude')
                    save_h5(hf, data=np.array([float(stLat)]), target='station_latitude')
                    save_h5(hf, data=np.array([float(stLon)]), target='station_longitude')
                    save_h5(hf, data=np.array([float(stElv)]), target='station_elevation')
                                
                    hf.close()
                    
                    del A, V 
            

            
            except:
                print('Unable to detect P arrival.')

      
   
    del tr
    
    
    
    
    