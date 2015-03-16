#!/usr/bin/python

import sys
import numpy as np
#mine
import utils

utils.usage(["<columns file>"],1,1)

colf = sys.argv[1]

if ".col" not in colf:
    print "run on a column file from outcar2columns.py"
    exit(0)

def grabV(colf):
    d = open(colf)
    [d.readline() for i in range(3)]
    l = d.readline()
    V = float(l.split()[-1])
    n = int(l.split()[0].split("=")[-1].strip(","))
    return V*n

def grabT(colf):
    d = open(colf)
    [d.readline() for i in range(8)]
    return float(d.readline().split()[1])

def parseCol(cold):
    return zip(*[map(float,x.split()[8:14]) for x in cold[7:] if len(x.split())==14 and x[0]!="S"])

V = grabV(colf)
T = grabT(colf)
coef = V/(5.*T*8.6173324E-5)
rho = 1/V
pxx,pyy,pzz,pxy,pxz,pyz = parseCol(open(colf).readlines())

pxx,pyy,pzz = np.array(pxx),np.array(pyy),np.array(pzz)
pxx = 0.5*(pxx-pzz)
pyy = 0.5*(pyy-pzz)
N=len(pxy)

mx=N/50
ts=range(1,mx)
tt=[]
for t in ts:
    tt.append(coef*sum([pxy[i]*pxy[i+t] + pxz[i]*pxz[i+t] + pyz[i]*pyz[i+t] + pxx[i]*pxx[i+t]+ pyy[i]*pyy[i+t] for i in range(0,N-mx+1,mx)]))
import pylab as pl
ts=[i*5./1000 for i in ts]
pl.plot(ts,tt/tt[0])
pl.show()
