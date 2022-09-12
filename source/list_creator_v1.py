# -*- coding: utf-8 -*-
"""
Created on Sun Aug  7 19:00:00 2022
Generates Compiled_list_v1.dat

@author: Alexander Nikolaev

"""

from datetime import datetime
#import re

# --- Function reading lists ---

def io(path):
    file = open(path,'r')
    data = file.readlines()
    file.close()
    return data

# --- Reading lists ---

path = ['f:\Injections_vs_DIP\lists\Motoba_2021.dat',
        'f:\Injections_vs_DIP\lists\Ohtani_2018.dat',
        'f:\Injections_vs_DIP\lists\Ohtani_2020_sharp.dat',
        'f:\Injections_vs_DIP\lists\Ohtani_2020_grad.dat']

#   --- Additional variables ---

head = ["%", "#"]
dt = []
dtSH = []
dtGR = []
info = []
infoSH = []
infoGR = []
output = []
outputGR = []
outputSH = []
index = []
check = []

# --- Motoba 2021 list ---

data = io(path[0])
#   SC YYYY MM DD HH MI SE  Lval  Rxy  MLT  MLat  Group
# ['$]YYYY MM DD HH MN]   [UNIXtime]   [VAPID  VAPMLT  RMLat]   [GID  GMLT  GMlat]   [E/P]   [gr/sh]'

for i in range(len(data)):
    try:
        data[i].index(head[1])
    except ValueError:
        [rbsp, datet, lval, rxy, rmlt, rmlat, ep] = data[i].split()
        dtt = datetime.strptime(datet, '%Y-%m-%d/%H:%M:%S')
        dtt = dtt.strftime('%Y %m %d %H %M')
        dt.append(datetime.strptime(dtt, '%Y %m %d %H %M'))
        info.append(rbsp + str('{0:7.2f}'.format(float(rmlt))) + '   ' + rmlat + '   ' + '#G' + '   ' + 'NaN' + '   ' + 'NaN' + '   ' + ep + '   ' + ' NaN')
    else:
         continue

# --- Ohtani [2018] list ---

# RBSP-A/B Year Month DayofMonth Hr Mn X_gsm Y_gsm Z_gsm MLat MLT L

data = io(path[1])

for i in range(len(data)-1):
    try:
        data[i].index(head[0])
    except ValueError:
        [rbsp, yy, mm, dd, hh, mn, x, y, z, mlat, mlt, L] = data[i].split()
        dt.append(datetime.strptime(yy+'/'+mm+'/'+dd+'T'+hh+':'+mn, '%Y/%m/%dT%H:%M'))
        info.append(rbsp + str('{0:7.2f}'.format(float(mlt))) + '   ' + mlat + '   ' + '#G' + '   ' + 'NaN' + '   ' + 'NaN' + '   ' + ' NaN' + '   ' + ' NaN')
    else:
        continue

# --- Ohtani 2020 list sharp---

data = io(path[2])

for i in range(len(data)-1):
    try:
        data[i].index(head[1])
    except ValueError:
        [date, time, goes, gmlat, gmlt, rbsp, rmlat, rmlt, dummy1, dummy2] = data[i].split()
        
        dt.append(datetime.strptime(date+'T'+time, '%Y/%m/%dT%H:%M'))
        dtSH.append(datetime.strptime(date+'T'+time, '%Y/%m/%dT%H:%M'))
        info.append(rbsp + str('{0:7.2f}'.format(float(rmlt))) + '   ' + rmlat + '   ' + goes + str('{0:6.2f}'.format(float(gmlt))) + '   ' + gmlat + '   ' + ' NaN' + '   ' + 'sharp')
        infoSH.append(rbsp + str('{0:7.2f}'.format(float(rmlt))) + '   ' + rmlat + '   ' + goes + str('{0:6.2f}'.format(float(gmlt))) + '   ' + gmlat + '   ' + ' NaN' + '   ' + 'sharp')
        
    else:
        continue

# --- Ohtani 2020 list gradual---

data = io(path[3])

for i in range(len(data)-1):
    try:
        data[i].index(head[1])
    except ValueError:
        [date, time, goes, gmlat, gmlt, rbsp, rmlat, rmlt, dummy1, dummy2] = data[i].split()
        
        dt.append(datetime.strptime(date+'T'+time, '%Y/%m/%dT%H:%M'))
        dtGR.append(datetime.strptime(date+'T'+time, '%Y/%m/%dT%H:%M'))
        info.append(rbsp + str('{0:7.2f}'.format(float(rmlt))) + '   ' + rmlat + '   ' + goes + str('{0:6.2f}'.format(float(gmlt))) + '   ' + gmlat + '   ' + ' NaN' + '   ' + 'gradD')
        infoGR.append(rbsp + str('{0:7.2f}'.format(float(rmlt))) + '   ' + rmlat + '   ' + goes + str('{0:6.2f}'.format(float(gmlt))) + '   ' + gmlat + '   ' + ' NaN' + '   ' + 'gradD')
        
    else:
        continue

#   --- Compile list ---

for i in range(len(dt)):
    if dt[i].strftime('%Y %m %d %H %M') not in check:
        check.append(dt[i].strftime('%Y %m %d %H %M'))
        output.append(dt[i].strftime('%Y %m %d %H %M')+'   '+str(dt[i].timestamp())+'   '+info[i])

#   --- Separate lists for sharp and gradual events ---

for i in range(len(dtGR)):
    outputGR.append(dtGR[i].strftime('%Y %m %d %H %M')+'   '+str(dtGR[i].timestamp())+'   '+infoGR[i])

for i in range(len(dtSH)):
    outputSH.append(dtSH[i].strftime('%Y %m %d %H %M')+'   '+str(dtSH[i].timestamp())+'   '+infoSH[i])        

#   --- Print output ---
output.sort()

file = open('f:\Injections_vs_DIP\lists\Compiled_list_v1.dat','w')
file.write('YYYY MM DD HH MN    UNIXtime   VAPID VAPMLT RMlat  GID  GMLT  GMlat   E/P   gr/sh')
file.write('\n')
for line in output:
    file.write(line)
    file.write('\n')
file.close()

"""
file = open('f:\Injections_vs_DIP\lists\Compiled_GRlist_v1.dat','w')
file.write('YYYY MM DD HH MN    UNIXtime   VAPID VAPMLT RMlat  GID  GMLT  GMlat   E/P   gr/sh')
file.write('\n')
for line in outputGR:
    file.write(line)
    file.write('\n')
file.close()

file = open('f:\Injections_vs_DIP\lists\Compiled_SHlist_v1.dat','w')
file.write('YYYY MM DD HH MN    UNIXtime   VAPID VAPMLT RMlat  GID  GMLT  GMlat   E/P   gr/sh')
file.write('\n')
for line in outputSH:
    file.write(line)
    file.write('\n')
file.close()
"""