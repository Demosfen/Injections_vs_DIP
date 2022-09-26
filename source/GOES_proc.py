# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 10:09:15 2022
Calculates S/C data for events lists

@author: Alexander
"""

import os
from datetime import datetime, timedelta
from geopack import geopack, t89
import matplotlib.pyplot as plt
import numpy as np
import pathlib
import urllib.request
import requests
import netCDF4

plt.style.use('seaborn-whitegrid')
#from math import sqrt

# ========================== Defs =================================
    
def file_exists(location):
    request = requests.head(location)
    return request.status_code == requests.codes.ok

# ------------------------------------------
"""
def cdfDownload(url):
    if file_exists(url):
        print("Downloading...")
        urllib.request.urlretrieve(url, goesFilePathLoc)
    return print("Download complete! Created file: "+goesFilePathLoc)
 """
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

# ================ PATH TO EVENTS LIST ============================

path = [pathCurrentFolder+'output/Events_MPB_GRlist_v1.dat',
        pathCurrentFolder+'data/GOES/']

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
if not os.path.exists(pathCurrentFolder+'data/GOES/variations'): 
    print("Create data/GOES/variations folder...")
    os.makedirs(pathCurrentFolder+'data/GOES/variations')

# --------- Define some variables -----------

head = ["YYYY"]
dt = []
dtime = []
output = []
mpb =[]
XGSM_min = [] ; YGSM_min = [] ; ZGSM_min = []
HP = [] ; HE = [] ; HN = []
XGSE = []; YGSE = []; ZGSE = []
time_field = []
HXGOES = [] ; HYGOES = [] ; HZGOES = []
time = []
plotTime = list(range(-600, 1800))
baseTime = datetime(2000,1,1,12,0)

# -------------- Reading events list --------------------

events = io(path[0])

for i in range(len(events)):
    try:
        events[i].index(head[0])
    except ValueError:
        [YYYY, MM, DD, HH, MN, dt_event, VAPID, VAPMLT, RMLat, GID, GMLT, GMLat, ep, grsh, dtmpb, dmpb] = events[i].split()
        
        # --------------- Paths ------------------
        
        goesYearPath = path[1]+YYYY
        goesMonthPath = goesYearPath + '/' + MM
        goesFilePath = goesMonthPath + '/' +'g'+GID.lower()+YYYY+'_magneto_512ms_'+YYYY+MM+DD+'_'+YYYY+MM+DD+'.csv'
        goesFilePathLoc = goesMonthPath + '/dn_goes-l2-orb1m_g'+GID+'_y'+YYYY+'_v0_0.nc'
        url = 'https://satdat.ngdc.noaa.gov/sem/goes/data/full/'+YYYY+'/'+MM+'/'+'goes'+GID.lower()+'/csv/g'+GID.lower()+'_magneto_512ms_'+YYYY+MM+DD+'_'+YYYY+MM+DD+'.csv'
        url_loc = 'https://satdat.ngdc.noaa.gov/sem/goes/data/sat_locations/goes'+GID+'/dn_goes-l2-orb1m_g'+GID+'_y'+YYYY+'_v0_0.nc'
       
        # ---------- Checking/Downloading RBSP file ------------
        
        if not os.path.exists(goesYearPath):
            print("*** "+YYYY+'/'+MM+'/'+DD+" ***")
            print("Creating year folder...")
            os.makedirs(path[1]+YYYY)
            print("Creating month folder...")
            os.makedirs(goesMonthPath)
            urllib.request.urlretrieve(url, goesFilePath)
            urllib.request.urlretrieve(url_loc, goesFilePathLoc)
        elif not os.path.exists(goesMonthPath):
            print("*** "+YYYY+'/'+MM+'/'+DD+" ***")
            print("Creating month folder...")
            os.makedirs(goesMonthPath)
            urllib.request.urlretrieve(url, goesFilePath)
            urllib.request.urlretrieve(url_loc, goesFilePathLoc)
        else:
            print("*** "+YYYY+'/'+MM+'/'+DD+" ***")
            print("Data folders for event "+YYYY+'/'+MM+'/'+DD+" found...")
            print("Searching data file...")
        
        if not os.path.exists(goesFilePath):
            print("GOES data file doesn't exist...")
            print("Searching...")
            urllib.request.urlretrieve(url, goesFilePath)
        else:
            print("Data file for "+YYYY+"/"+MM+"/"+DD+" found...")
            
        if not os.path.exists(goesFilePathLoc):
            print("GOES location file doesn't exist...")
            print("Searching...")
            urllib.request.urlretrieve(url_loc, goesFilePathLoc)
        else:
            print("Location file for "+YYYY+"/"+MM+"/"+DD+" found...\n")
            
# -------------- Reading data --------------------            
            
        goes_field = io(goesFilePath)
        locationDataset = netCDF4.Dataset(goesFilePathLoc)
        xyzGSENC = locationDataset['gse_xyz']
        timeLocNC = locationDataset['time']
        timeLoc = timeLocNC[:]
        locationTime = [None] * len(timeLoc)
        xyzGSE = xyzGSENC[:]/6371.15
        for k in range(len(timeLoc)): locationTime[k] = (baseTime+timedelta(seconds=timeLoc[k])).timestamp()

# -------------- Reading events list --------------------

        for j in range(596,len(goes_field)):
            k = j-596
            [DTime, 
             BX1Q,BX1,BY1Q,BY1,BZ1Q,BZ1,
             BXSC1Q,BXSC1,BYSC1Q,BYSC1,BZSC1Q,BZSC1,BTSC1Q,BTSC1,
             BX2Q,BX2,BY2Q,BY2,BZ2Q,BZ2,
             BXSC2Q,BXSC2,BYSC2Q,BYSC2,BZSC2Q,BZSC2,BTSC2Q,BTSC2,
             HP1Q,HP1,HE1Q,HE1,HN1Q,HN1,HT1Q,HT1,
             HP2Q,HP2,HE2Q,HE2,HN2Q,HN2,HT2Q,HT_2] = goes_field[j].split(',')
            
            HP.append(float(HP2))
            HE.append(float(HE2))
            HN.append(float(HN2))
            time_field.append(datetime.strptime(DTime, '%Y-%m-%d %H:%M:%S.%f').timestamp())
            continue
        
        dtime0 = datetime.strptime(YYYY+'-'+MM+'-'+DD+" 00:00:00", '%Y-%m-%d %H:%M:%S').timestamp()
        i_t0 = locationTime.index(dtime0)
        dtime1 = datetime.strptime(YYYY+'-'+MM+'-'+DD+" 23:59:00", '%Y-%m-%d %H:%M:%S').timestamp()
        x = np.arange(dtime0,dtime1,60,'float')
        
        HP_min = np.interp(x, time_field,HP)
        HE_min = np.interp(x, time_field,HE)
        HN_min = np.interp(x, time_field,HN)
        
        for j in range(len(HP_min)):
            ps = geopack.recalc(x[j])
            XGSM, YGSM, ZGSM = geopack.gsmgse(xyzGSE[i_t0+j,0], xyzGSE[i_t0+j,1], xyzGSE[i_t0+j,2], -1)
            bxigrf,byigrf,bzigrf = geopack.igrf_gsm(XGSM,YGSM,ZGSM)
            heigrf,hpigrf,hnigrf = geopack.bcarsp(XGSM,YGSM,ZGSM,bxigrf,byigrf,bzigrf)
            bxt89,byt89,bzt89 = t89.t89(2,ps,XGSM,YGSM,ZGSM)
            het89,hpt89,hnt89 = geopack.bcarsp(XGSM,YGSM,ZGSM,bxt89,byt89,bzt89)
            r,theta,phi = geopack.sphcar(XGSM,YGSM,ZGSM,-1)
            BXGEO,BYGEO,BZGEO = geopack.bspcar(theta, phi, HE_min[j], HP_min[j]*(-1), HN_min[j])
            BXGSM,BYGSM,BZGSM = geopack.geogsm(BXGEO, BYGEO, BZGEO, 1)
            HXGSM,HZGSM,HYGSM = geopack.bcarsp(XGSM,YGSM,ZGSM, BXGSM,BYGSM,BZGSM)
            HXGOES.append(HXGSM-(het89+heigrf))
            HYGOES.append(HYGSM-(hnt89+hnigrf))
            HZGOES.append(HZGSM*(-1)-(hpt89*(-1)+hpigrf*(-1)))
            XGSM_min.append(XGSM) ; YGSM_min.append(YGSM) ; ZGSM_min.append(ZGSM)

        for j in range(len(HP_min)): output.append(datetime.fromtimestamp(x[j]).strftime('%Y-%m-%d/%H:%M:%S')+' '+str(HXGOES[j])+' '+str(HYGOES[j])+' '+str(HZGOES[j]))
        
        fig = plt.figure()
        ax = plt.axes()
        ax.plot(x,HZGOES)
        plt.show()
        
        file = open('f:/#Research/Injections_vs_DIP/data/GOES/variations/var_'+YYYY+MM+DD+'.dat','w')
        file.write('Epoch               HXGSM  HYGSM  HZGSM')
        file.write('\n')
        for line in output:
            file.write(line)
            file.write('\n')
        file.close()
        
        del HP[:], HE[:], HN[:], HXGOES[:], HYGOES[:], HZGOES[:]
        del XGSM_min[:], YGSM_min[:], ZGSM_min[:]
        del time_field[:]
        del output[:]
        
    else:
        continue
    