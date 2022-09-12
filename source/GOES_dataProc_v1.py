# -*- coding: utf-8 -*-
"""
Created on Sat Sep  3 18:12:03 2022

@author: Alexander
"""

from datetime import datetime
from geopack import geopack, t89
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
plt.style.use('seaborn-whitegrid')
#from math import sqrt

# --- Function reading lists ---

def io(path):
    file = open(path,'r')
    data = file.readlines()
    file.close()
    return data

# --- Reading lists ---

def savitzky_golay(y, window_size, order, deriv=0, rate=1):

    import numpy as np
    from math import factorial
    
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

path = ['f:/#Research/Injections_vs_DIP/output/Events_MPB_GRlist_v1.dat',
        'f:/#Research/Injections_vs_DIP/data/GOES/']
head = ["YYYY"]
dt = []
dtime = []
time_field = []
time_loc = []
output = []
XGSM = []
YGSM = []
ZGSM = []
HP = []
HE = []
HN = []
HXGOES = []
HYGOES = []
HZGOES = []
plotTime = list(range(-600, 1800))

# --- Reading events list ---

events = io(path[0])

for i in range(len(events)):
    try:
        events[i].index(head[0])
    except ValueError:
        [YYYY, MM, DD, HH, MN, dt_event, VAPID, VAPMLT, RMLat, GID, GMLT, GMLat, ep, grsh, dtmpb, dmpb] = events[i].split()

        # --- Reading GOES magnetic file ---

        goes_field = io(path[1]+YYYY+MM+DD+'_g'+GID+'.csv')
        goes_loc = io(path[1]+YYYY+MM+DD+'_g'+GID+'_loc.dat')
        
        for j in range(596,len(goes_field)):
            k = j-596
            [DTime, 
             BX1Q,BX1,BY1Q,BY1,BZ1Q,BZ1,
             BXSC1Q,BXSC1,BYSC1Q,BYSC1,BZSC1Q,BZSC1,BTSC1Q,BTSC1,
             BX2Q,BX2,BY2Q,BY2,BZ2Q,BZ2,
             BXSC2Q,BXSC2,BYSC2Q,BYSC2,BZSC2Q,BZSC2,BTSC2Q,BTSC2,
             HP1Q,HP1,HE1Q,HE1,HN1Q,HN1,HT1Q,HT1,
             HP2Q,HP2,HE2Q,HE2,HN2Q,HN2,HT2Q,HT_2] = goes_field[j].split(',')
            
            #DTime=DTime[:19]
            #dtime.append(DTime)
            HP.append(float(HP2))
            HE.append(float(HE2))
            HN.append(float(HN2))      
            time_field.append(datetime.strptime(DTime, '%Y-%m-%d %H:%M:%S.%f').timestamp())
            #if time[k] == float(dt_event): i_t0 = k                
            
            continue
        dtime0 = datetime.strptime(YYYY+'-'+MM+'-'+DD+" 00:00:00", '%Y-%m-%d %H:%M:%S').timestamp()
        dtime1 = datetime.strptime(YYYY+'-'+MM+'-'+DD+" 23:59:00", '%Y-%m-%d %H:%M:%S').timestamp()
        x = np.arange(dtime0,dtime1,60,'float')
        
        HP_min = np.interp(x, time_field,HP)
        HE_min = np.interp(x, time_field,HE)
        HN_min = np.interp(x, time_field,HN)
        
        for j in range(1, len(goes_loc)):
            [Date, Time, X, Y, Z] = goes_loc[j].split()
            DTime=Date+'T'+Time
            time_loc.append(datetime.strptime(DTime, '%d-%m-%YT%H:%M:%S.%f').timestamp())
            XGSM.append(float(X))
            YGSM.append(float(Y))
            ZGSM.append(float(Z))
            continue
        XGSM_min = np.interp(x, time_loc,XGSM)
        YGSM_min = np.interp(x, time_loc,YGSM)
        ZGSM_min = np.interp(x, time_loc,ZGSM)
        
        #HP_sm=savitzky_golay(HP_min, 61, 3)
        #HE_sm=savitzky_golay(HE_min, 61, 3)
        #HN_sm=savitzky_golay(HN_min, 61, 3)
        
        for j in range(len(HP_min)):
            ps = geopack.recalc(x[j])
            bxigrf,byigrf,bzigrf = geopack.igrf_gsm(XGSM_min[j],YGSM_min[j],ZGSM_min[j])
            heigrf,hpigrf,hnigrf = geopack.bcarsp(XGSM_min[j],YGSM_min[j],ZGSM_min[j],bxigrf,byigrf,bzigrf)
            bxt89,byt89,bzt89 = t89.t89(2,ps,XGSM_min[j],YGSM_min[j],ZGSM_min[j])
            het89,hpt89,hnt89 = geopack.bcarsp(XGSM_min[j],YGSM_min[j],ZGSM_min[j],bxt89,byt89,bzt89)
            r,theta,phi = geopack.sphcar(XGSM_min[j],YGSM_min[j],ZGSM_min[j],-1)
            BXGEO,BYGEO,BZGEO = geopack.bspcar(theta, phi, HE_min[j], HP_min[j]*(-1), HN_min[j])
            BXGSM,BYGSM,BZGSM = geopack.geogsm(BXGEO, BYGEO, BZGEO, 1)
            HXGSM,HZGSM,HYGSM = geopack.bcarsp(XGSM_min[j],YGSM_min[j],ZGSM_min[j], BXGSM,BYGSM,BZGSM)
            #HXGOES.append(HXGSM-(het89+heigrf))
            #HYGOES.append(HYGSM-(hnt89+hnigrf))
            HZGOES.append(HZGSM*(-1)-(hpt89*(-1)+hpigrf*(-1)))
            #HXGOES.append(HXGSM)
            #HYGOES.append(HYGSM)
            #HZGOES.append(HP_min[j])
            
            continue
        
        for j in range(len(HP_min)): output.append(datetime.fromtimestamp(x[j]).strftime('%Y-%m-%dT%H:%M:%S')+' '+str(HZGOES[j]))
        
        fig = plt.figure()
        ax = plt.axes()
        ax.plot(x,HZGOES)
        plt.show()
        
        file = open('f:/#Research/Injections_vs_DIP/data/GOES/var_'+YYYY+MM+DD+'.dat','w')
        #file = open('f:/#Research/Injections_vs_DIP/data/GOES/var_'+YYYY+MM+DD+'_noT89.dat','w')
        file.write('Epoch               HZGSM')
        file.write('\n')
        for line in output:
            file.write(line)
            file.write('\n')
        file.close()
        
        del HP[:], HE[:], HN[:], HXGOES[:], HYGOES[:], HZGOES[:]
        del XGSM[:], YGSM[:], ZGSM[:]
        del time_field[:], time_loc[:]
        del output[:]

    else:
        continue
    