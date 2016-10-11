#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  ion_motion_track.py
#  
#  Copyright 2016 weiwei <weiwei@xps8700>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  



import sys
import math
import pylab
import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as pd
import heapq

box = [0,0,0]
ionid = []
nid = []
m = 0
si = 0
ins = []
insl = []
ionpn = []
ionpo = []
dx = []
dy = []
dz = []
ddt = []
nsmallestt = 1
ddsm = []
nstep = 0

pos = pd.DataFrame()

# range of solvation shell
ss = 6.1


def toprime(pos,box):
    posn = []
    for i in range(3):
        n1 = math.floor((pos[i]-box[i*2])/(box[i*2+1]-box[i*2]))
        n2 = math.floor((pos[i]-box[i*2+1])/(box[i*2+1]-box[i*2]))
        n3 = n1 if n1 > n2 else n2
        posn.append(pos[i] - n3*(box[i*2+1]-box[i*2]))
    for i in range(3):
        if posn[i] > box[2*i+1] or posn[i] < box[2*i]:
            print "move to prime is not successful"
    return posn

def dist(x,y):
    dx = x[0]-y[0]
    dy = x[1]-y[1]
    dz = x[2]-y[2]
    return math.sqrt(dx*dx+dy*dy+dz*dz)
    
f = open("npt.lammpstrj",'r')
lines = f.readlines()
for i in range(len(lines)):
    li = lines[i].split()
    if (len(li) == 2 and li[1] == "TIMESTEP"):
        m = m+1
        if m == 2:
            break
    if (len(li) == 6 and li[0].isdigit()):
        if int(li[2]) == 18:
            ionid.append(int(li[0]))
        if int(li[2]) == 15:
            nid.append(int(li[0]))

ionid.sort()
nid.sort()
if len(nid) != len(ionid):
    print "Number of ions is not the same as nitrogen"

lxx = np.arange(0,len(nid),1)
dicn = dict(zip(nid,lxx))

si = random.randint(ionid[0],ionid[-1])
print si
#nid.append(si)         
for i in range(len(lines)):
    li = lines[i].split()
    if (len(li)==6 and li[1] == "BOX"):
        nstep=nstep+1
        boxt = [float(lines[i+1].split()[0]),float(lines[i+1].split()[1]),float(lines[i+2].split()[0]),float(lines[i+2].split()[1]),float(lines[i+3].split()[0]),float(lines[i+3].split()[1])]
        box[0] = float(lines[i+1].split()[1]) - float(lines[i+1].split()[0])
        box[1] = float(lines[i+2].split()[1]) - float(lines[i+2].split()[0])
        box[2] = float(lines[i+3].split()[1]) - float(lines[i+3].split()[0])
        if len(pos) != 0:
            for k in range(len(pos)):
                pn = [pos.ix[k,'x'],pos.ix[k,'y'],pos.ix[k,'z']]
                pn = toprime(pn,boxt)
                pi = toprime(ionpn,boxt)
                d = dist(pn,pi)
                #print pos.ix[k,'type'],pos.ix[k,'x'],pos.ix[k,'y'],pos.ix[k,'z']
                #d = math.sqrt((pos.ix[k,'x']-ionpn[0])*(pos.ix[k,'x']-ionpn[0])+(pos.ix[k,'y']-ionpn[1])*(pos.ix[k,'y']-ionpn[1])+(pos.ix[k,'z']-ionpn[2])*(pos.ix[k,'z']-ionpn[2]))
               # print d                
                ddt.append(d)
            #get the smallest three element in ddt,it is in dds
            dds = heapq.nsmallest(nsmallestt,ddt)
            ddsm.append(np.mean(np.asarray(dds)))
            for kk in range(nsmallestt):
                ins.append(int(dicn[int(pos.ix[ddt.index(dds[kk]),'atomID'])]))               
                #if d <= ss:                 
                    ##print "in solvaton", pos.ix[k,'atomID']
                    #ins.append(int(dicn[int(pos.ix[k,'atomID'])]))
            insl.append(np.asarray(ins))
            print nstep,ins,ddsm[-1]
            ins = [] 
            dds = []
            ddt = []                   
            pos = pd.DataFrame()
            
    #1. track the motion of ions
    #2. find the list of nitrogens within ion solvation shell
    if (len(li)==6 and li[0].isdigit()):
        if int(li[0]) == si:
            ionpn=[float(li[3]),float(li[4]),float(li[5])]
            if len(ionpo) == 0:
                ionpo = ionpn
            else:
                dx.append(ionpn[0]-ionpo[0])
                dy.append(ionpn[1]-ionpo[1])
                dz.append(ionpn[2]-ionpo[2])
            #print ionpn
        if int(li[2]) == 15:
            li = map(lambda x:float(x),li)
            temp = pd.Series(li,index = ['atomID','chainID','type','x','y','z'])
            #print temp.values
            pos = pos.append([temp],ignore_index=True)
        
x = np.arange(0,len(dx),1)
plt.subplot(3,1,1)
plt.plot(x,dx,label='dx')
plt.plot(x,dy,label='dy')
plt.plot(x,dz,label='dz') 
plt.legend()  
plt.subplot(3,1,2)
plt.plot(x,np.asarray(insl)[:,0],label='coord')
#plt.xlabel('Coordination VS time')
plt.legend() 
plt.subplot(3,1,3)
plt.plot(x,ddsm,label='dist')
plt.xlabel('t/10ps')
plt.legend() 
plt.show()
    
        
    

