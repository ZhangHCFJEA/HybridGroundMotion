# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 09:46:31 2021

@author: ZhangHC
"""
import os
import h5py
import tensorflow as tf
from sklearn.model_selection import train_test_split
import numpy as np
from keras.models import load_model
import matplotlib.pyplot as plt
from obspy.core import read, UTCDateTime, Stream
from obspy.geodetics.base import gps2dist_azimuth as g2d

path = 'F:/kiknet_ml'
fileName = 'kiknet.hdf5'
f = h5py.File(os.path.join(path, fileName), 'r')


loc, ID = [], []
n = 0
while n < len(f['recordID']):
    print(n)
    flag = f['recordID'][n]
    ID.append(str(f['recordID'][n][6:16], encoding = "utf8"))
    loc.append([f['event_latitude'][n], f['event_longitude'][n], f['focal depth'][n], f['Mag'][n]])
    n += 1

np.savetxt('events.txt', loc)    
filename = open('id.txt', 'w')
for value in ID:
    filename.write(value)
    filename.write('\n')
filename.close()    



count, loc = [], []
flag = '0503201053'
for n in range(len(f['recordID'])):
    if str(f['recordID'][n][6:16], encoding = "utf8") == flag:
        count.append(n)
        loc.append([f['station_latitude'][n], f['station_longitude'][n]])
        print(f['event_latitude'][n])
        print(f['event_longitude'][n])
        print(f['recordID'][n])

# for i in range(len(loc)):
#     plt.plot(loc[i][1], loc[i][0], 'b^')
# plt.plot(130.175, 33.738, 'r+')

# originTime = UTCDateTime('2005-03-20T10:53:40.32')
originTime = UTCDateTime('2005-03-20T01:53:40.32')
epc = [33.731967, 130.176333, 9.24]
recordPath = 'F:/kiknet_ml/examples/event1'
records = os.listdir(recordPath)
P = np.loadtxt('F:/kiknet_ml/examples/Ppick1.txt')
tr = Stream()
for r in range(2, len(records), 3):
    tr += read(os.path.join(recordPath, records[r]))
    index = r//3
    stDist = g2d(epc[0], epc[1], tr[index].stats.knet.stla, tr[index].stats.knet.stlo)
    tr[index].data = tr[index].data-np.mean(tr[index].data[:1000])
    tr[index].stats.distance = stDist[0]
    # tr[index].data = tr[index].data/max([np.abs(tr[index].data.min()), tr[index].data.max()])+stDist[0]/1000.

tr.trim(starttime=originTime, endtime=originTime+40 )
fig = plt.figure()
tr.plot(type='section', plot_dx=20e3, recordlength=40,
        time_down=True, linewidth=.25, grid_linewidth=.25, show=False, fig=fig)


#
model1 = '1best_model.54-0.41.h5'
model2 = '2best_model.43-0.31.h5'
model3 = '3best_model.56-0.27.h5'
model4 = '4best_model.59-0.24.h5'
model5 = '5best_model.56-0.22.h5'
model6 = '6best_model.59-0.21.h5'
model7 = '7best_model.43-0.19.h5'
model8 = '8best_model.56-0.18.h5'
model9 = '9best_model.34-0.18.h5'
model10 = '10best_model.59-0.16.h5'
m1 = load_model(os.path.join(path, model1))
m2 = load_model(os.path.join(path, model2))
m3 = load_model(os.path.join(path, model3))
m4 = load_model(os.path.join(path, model4))
m5 = load_model(os.path.join(path, model5))
m6 = load_model(os.path.join(path, model6))
m7 = load_model(os.path.join(path, model7))
m8 = load_model(os.path.join(path, model8))
m9 = load_model(os.path.join(path, model9))
m10 = load_model(os.path.join(path, model10))

#
X = f['accData']
Y = f['JMA intensity']
tf.random.set_seed(1)
sampling_rate = 100.
# train_x, test_x, train_Y, test_Y = train_test_split(np.array(X[:, :, :]), np.array(Y), test_size=0.25, random_state=7)
flag = '0503201053'
for win in range(10):
    print(win)
    waveLen = int((win+1)*sampling_rate)
    train_X = X[:, :, :waveLen]
    train_X = np.array(train_X).transpose(0, 2, 1)
    train_X = train_X/np.std(train_X)
    
    #
    if win == 0:
        m = m1
    elif win == 1:
        m = m2
    elif win == 2:
        m = m3
    elif win == 3:
        m = m4
    elif win == 4:
        m = m5
    elif win == 5:
        m = m6
    elif win == 6:
        m = m7
    elif win == 7:
        m = m8
    elif win == 8:
        m = m9
    else:
        m = m10
    pred_Y = m.predict(train_X)
    
    #
    output, stList = [], []
    for n in range(len(f['recordID'])):
        if str(f['recordID'][n][6:16], encoding = "utf8") == flag:
            output.append([pred_Y[n][0], Y[n]])
            stList.append([f['station_latitude'][n], f['station_longitude'][n]])
            # print([f['station_latitude'][n], f['station_longitude'][n]])
    
    np.savetxt(str(win)+'.txt', output, fmt='%f', delimiter=' ' )
    del train_X, pred_Y, output

np.savetxt('stList.txt', stList, fmt='%f', delimiter=' ' )















