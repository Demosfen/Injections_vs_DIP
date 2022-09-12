#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 18:21:15 2022

@author: alexander
"""

from datetime import datetime, timedelta
from math import sqrt

# --- Function reading lists ---
def io(path):
    file = open(path,'r')
    data = file.readlines()
    file.close()
    return data

# --- Reading Events_v0 ---
path = ['/media/alexander/Botanica/#Research/Injections_vs_DIP/lists/Events_v0.dat']
head = "$"
dt = []
info = []
output = []

events = io(path[1])
mpb = io(path[0])

# --- MPB onset list ---
j = 0
j1 = 0
for i in range(len(events)-1):
    try:
        events[i].index(head)
    except ValueError: