# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 16:03:46 2021

@author: ZhangHC
"""
import os
import numpy as np

def JMATravelTime(model, depth, epiDist):
    '''
    ***** Travel time table downloaded from JMA website. *****
    Table look up to provide seismic ray travel times in Japan.
    NOTE: JUST PROVIDE THE CLOSEST VALUES IN THE TABLE, NO INTERPOLATION!
    
    ----------
    Useage:
        JMATravelTime(model, depth, epiDist)
    
    Parameters
    ----------
    model : string
        JMA support model names.
        Include: tjma2001, tll2001, tt83a, ttjb, ttll
    depth : float
        Source depth, in km.
    epiDist : float
        Epicentral distance, in km.

    Returns
    -------
    travelTime : list
        Returns P and S wave travle times, in second.

    '''
    dataPath = 'D:/桌面/数据/JMA/arrivalTime'
    if model == 'tjma2001' or model == 'tll2001' or model == 'tt83a' or model == 'ttjb' or model == 'ttll' :
        # travelTime = np.genfromtxt(os.path.join(dataPath, model, model), usecols=[1, 3, 4, 5])
        PTime = np.genfromtxt(os.path.join(dataPath, model, model), usecols=[1])
        STime = np.genfromtxt(os.path.join(dataPath, model, model), usecols=[3])
        dep = np.genfromtxt(os.path.join(dataPath, model, model), usecols=[4])
        dist = np.genfromtxt(os.path.join(dataPath, model, model), usecols=[5])
        
        tmp = np.abs(depth-dep)
        k1 = np.where(tmp == tmp.min())[0][0]
        index = np.argwhere(tmp==min(tmp))
        tmp = np.abs(epiDist-dist[index])
        k2 = np.where(tmp == tmp.min())[0][0]
        
        travelTime = [PTime[k1+k2], STime[k1+k2]]

    else:
        print('Model Not Found!')
        travelTime = [np.nan, np.nan]
    
    return travelTime