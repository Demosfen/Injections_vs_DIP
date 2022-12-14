# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 10:09:15 2022
Calculates S/C data for events lists

@author: Alexander
"""

import os
# https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/ - CDF library download here
os.environ["CDF_LIB"] = 'C:/Program Files/CDF_Distribution/cdf38_1-dist/'
from datetime import datetime
from geopack import geopack, t89
import matplotlib.pyplot as plt
import numpy as np
import pathlib
import urllib.request
import requests
# To install spacepy in miniconda type: pip install spacepy
import spacepy.pycdf as cdf

plt.style.use('seaborn-whitegrid')
#from math import sqrt

# ========================== Defs =================================
    
def file_exists(location):
    request = requests.head(location)
    return request.status_code == requests.codes.ok

# ------------------------------------------

def cdfDownload(url,cdfVersions):
    for ii in range(len(cdfVersions)):
        if file_exists(url+cdfVersions[ii]+'.cdf'):
            print("Downloading...")
            urllib.request.urlretrieve(url+cdfVersions[ii]+'.cdf', rbspFilePath)
            break
        else:
            continue
    return print("Download complete! Created file: "+rbspFilePath)
 
# --- Function reading lists ---

def io(path):
    file = open(path,'r')
    data = file.readlines()
    file.close()
    return data

# --- Smooth algorithm ---

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

pathCurrentFolder = str(pathlib.Path(__file__).parent.resolve())
pathLen = len(pathCurrentFolder)
pathCurrentFolder = pathCurrentFolder[:pathLen-6].replace("\\", '/')
path = [pathCurrentFolder+'output/Events_v1.dat',
        pathCurrentFolder+'data/RBSP/']

# =================================================================

# --------- Generate project folders -----------

if not os.path.exists(pathCurrentFolder+'data/'):
    print("Create /data folder...")
    os.makedirs(pathCurrentFolder+'data')
    os.makedirs(pathCurrentFolder+'data/MPB')
    os.makedirs(pathCurrentFolder+'data/RBSP')
    os.makedirs(pathCurrentFolder+'data/RBSP/variations')
    os.makedirs(pathCurrentFolder+'data/GOES')
    os.makedirs(pathCurrentFolder+'data/GOES/variations')
if not os.path.exists(pathCurrentFolder+'output/'):
    print("Create /output folder...")
    os.makedirs(pathCurrentFolder+'output')
if not os.path.exists(pathCurrentFolder+'lists/'): 
    print("Create /lists folder...")
    os.makedirs(pathCurrentFolder+'lists')
if not os.path.exists(pathCurrentFolder+'graph/'): 
    print("Create /graph folder...")
    os.makedirs(pathCurrentFolder+'graph')
if not os.path.exists(pathCurrentFolder+'data/RBSP/variations'): 
    print("Create data/RBSP/variations folder...")
    os.makedirs(pathCurrentFolder+'data/RBSP/variations')

# --------- Define some variables -----------

head = ["YYYY"]
dt = []
dtime = []
output = []
mpb =[]
XGSM = [] ; YGSM = [] ; ZGSM = []
BXRBSP = [] ; BYRBSP = [] ; BZRBSP = []
time = []
plotTime = list(range(-600, 1800))
cdfVersions = ['1.3.3', '1.3.4', '1.3.5', '1.3.6', '1.6.1','1.6.2','1.6.3', '1.7.1', '1.7.2', '1.7.3']

# -------------- Reading events list --------------------

events = io(path[0])

for i in range(len(events)):
    try:
        events[i].index(head[0])
    except ValueError:
        [YYYY, MM, DD, HH, MN, dt_event, VAPID, VAPMLT, RMLat, GID, GMLT, GMLat, ep, grsh, dtmpb, dmpb] = events[i].split()
        
        # --------------- Paths ------------------
        
        rbspYearPath = path[1]+YYYY
        rbspMonthPath = rbspYearPath + '/' + MM
        rbspFilePath = rbspMonthPath + '/' +'rbsp-'+VAPID.lower()+'_magnetometer_4sec-gsm_emfisis-l3_'+YYYY+MM+DD+'.cdf'
        url = 'https://cdaweb.gsfc.nasa.gov/pub/data/rbsp/rbsp'+VAPID.lower()+'/l3/emfisis/magnetometer/4sec/gsm/'+YYYY+'/rbsp-'+VAPID.lower()+'_magnetometer_4sec-gsm_emfisis-l3_'+YYYY+MM+DD+'_v'
       
        # ---------- Checking/Downloading RBSP file ------------
        
        if not os.path.exists(rbspYearPath):
            print("*** "+YYYY+'/'+MM+'/'+DD+" ***")
            print("Creating year folder...")
            os.makedirs(path[1]+YYYY)
            print("Creating month folder...")
            os.makedirs(rbspMonthPath)
            cdfDownload(url,cdfVersions)
        elif not os.path.exists(rbspMonthPath):
            print("*** "+YYYY+'/'+MM+'/'+DD+" ***")
            print("Creating month folder...")
            os.makedirs(rbspMonthPath)
            cdfDownload(url,cdfVersions)
        else:
            print("*** "+YYYY+'/'+MM+'/'+DD+" ***")
            print("Data folders for event "+YYYY+'/'+MM+'/'+DD+" found...")
            print("Searching data file...")
        
        if not os.path.exists(rbspFilePath):
            print("RBSP file doesn't exist...")
            print("Searching...")
            cdfDownload(url,cdfVersions)
        else:
            print("Data file for "+YYYY+"/"+MM+"/"+DD+" found...\n")
            
        rbspCDF = cdf.CDF(rbspFilePath)
        #print(rbspCDF)
        Epoch = rbspCDF['Epoch'][...]
        dataMag = rbspCDF['Mag'][...]
        dataCoord = rbspCDF['coordinates'][...]
        
        for dt in Epoch: time.append(dt.replace(second=0, microsecond=0).timestamp())
        i_t0 = time.index(float(dt_event))
        del time[:]
        for dt in Epoch: time.append(dt.replace(microsecond=0).timestamp())
        time = time[i_t0-900:i_t0+900]
        XGSM = dataCoord[i_t0-900:i_t0+900,0]/6371.15
        YGSM = dataCoord[i_t0-900:i_t0+900,1]/6371.15
        ZGSM = dataCoord[i_t0-900:i_t0+900,2]/6371.15
        BX = dataMag[i_t0-900:i_t0+900,0]
        BY = dataMag[i_t0-900:i_t0+900,1]
        BZ = dataMag[i_t0-900:i_t0+900,2]
        
        
        for j in range(len(XGSM)):
            ps = geopack.recalc(time[j])
            bxigrf,byigrf,bzigrf = geopack.igrf_gsm(XGSM[j],YGSM[j],ZGSM[j])
            bxigrf,bzigrf,byigrf = geopack.bcarsp(XGSM[j],YGSM[j],ZGSM[j],bxigrf,byigrf,bzigrf)
            bxt89,byt89,bzt89 = t89.t89(2,ps,XGSM[j],YGSM[j],ZGSM[j])
            bxt89,bzt89,byt89 = geopack.bcarsp(XGSM[j],YGSM[j],ZGSM[j],bxt89,byt89,bzt89)
            HX,HZ,HY = geopack.bcarsp(XGSM[j],YGSM[j],ZGSM[j],float(BX[j]),float(BY[j]),float(BZ[j]))
            BXRBSP.append(HX-(bxigrf+bxt89))
            BYRBSP.append(HY-(byigrf+byt89))
            BZRBSP.append(HZ*(-1)-(bzigrf*(-1)+bzt89*(-1)))
            continue

        for l in range(len(BZRBSP)): 
            dateTime = datetime.fromtimestamp(time[l]).strftime('%Y-%m-%d/%H:%M:%S')
            output.append('%s %.2f %.2f %.2f %.2f %.2f %.2f' % (dateTime, XGSM[l], YGSM[l], ZGSM[l], BXRBSP[l], BYRBSP[l], BZRBSP[l]))
        
        fig = plt.figure()
        ax = plt.axes()
        ax.plot(time,BZRBSP)
        #ax.plot(time,BZRBSP_sm[i_t0-600:i_t0+1800])
        plt.show()
        
        file = open(pathCurrentFolder+'data/RBSP/variations/VAP'+VAPID+'_'+YYYY+MM+DD+HH+MN+'.dat','w')
        file.write('YYYY-MM-DD/HH:MN:SS     X   Y    Z    BX    BY   BZ')
        file.write('\n')
        for line in output:
            file.write(line)
            file.write('\n')
        file.close()
        
        #BZRBSP_sm = np.array(BZRBSP_sm).tolist()
        XGSM = XGSM.tolist()
        YGSM = YGSM.tolist()
        ZGSM = ZGSM.tolist()
        del BXRBSP[:], BYRBSP[:], BZRBSP[:] #, BZRBSP_sm[:]
        del XGSM[:], YGSM[:], ZGSM[:]
        del time[:] #, dtime[:]
        del output[:]

    else:
        continue
    