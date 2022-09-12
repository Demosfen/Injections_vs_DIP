# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 12:41:01 2022
Calculates MPB data for events lists

@author: Alexander Nikolaev


"""

from datetime import datetime
#from math import sqrt

# --- Function reading lists ---

def io(path):
    file = open(path,'r')
    data = file.readlines()
    file.close()
    return data

# --- Reading lists ---

path = ['f:/Injections_vs_DIP/lists/Compiled_GRlist_v1.dat',
        'f:/Injections_vs_DIP/data/MPB/']
head = ["YYYY"]
dt = []
output = []
mpb =[]

# --- Reading events list ---

events = io(path[0])

for i in range(len(events)):
    try:
        events[i].index(head[0])
    except ValueError:
        [YYYY, MM, DD, HH, MN, dt_event, VAPID, VAPMLT, RMLat, GID, GMLT, GMLat, ep, grsh] = events[i].split()

# --- Reading MPB file ---

        mpb_arr = io(path[1]+YYYY+'.dat')
        [MYYYY, MMM, MDD, MHH, MMN, MMPB] = mpb_arr[1].split()
        dt_mpb0 = datetime.strptime(MYYYY+'/'+MMM+'/'+MDD+'T'+MHH+':'+MMN, '%Y/%m/%dT%H:%M').timestamp()
        dt0 = int((float(dt_event)-dt_mpb0)/60)-30
        dt1 = int((float(dt_event)-dt_mpb0)/60)+30
        mpb_arr = mpb_arr[dt0:dt1]
        mpb =[]
        for j in range(len(mpb_arr)-1):
            [MYYYY, MMM, MDD, MHH, MMN, MMPB] = mpb_arr[j].split()
            mpb.append(float(MMPB))
        mpb0 = min(mpb)
        index_mpb0 = min(range(len(mpb)), key=mpb.__getitem__)
        mpb1 = max(mpb)
        index_mpb1 = max(range(len(mpb)), key=mpb.__getitem__)
        dmpb = mpb1-mpb0
        dtmpb = index_mpb1-index_mpb0
        output.append((events[i] + '   ' + str(dtmpb) + '   ' + str('{0:7.2f}'.format(dmpb))).replace('\n',''))
        output
    else:
        continue

file = open('f:\Injections_vs_DIP\output\Events_MPB_GRlist_v1.dat','w')
file.write('YYYY MM DD HH MN    UNIXtime   VAPID VAPMLT RMlat  GID  GMLT  GMlat   E/P   gr/sh   dtMPB   dMPB')
file.write('\n')
for line in output:
    file.write(line)
    file.write('\n')
file.close()