#!/usr/bin/python

import sys
#mine
import utils

utils.usage(["<columns file>"],1,1)

colf = sys.argv[1]

if ".col" not in colf:
    print "run on a column file from outcar2columns.py"
    exit(0)


def parseCol(cold):
    return zip(*[map(float,x.split()[8:14]) for x in cold[7:] if len(x.split())==14 and x[0]!="S"])

pxx,pyy,pzz,pxy,pxz,pyz = parseCol(open(colf).readlines())

N=len(pxy)

ts=range(1,N/10)
tt=[]
for t in ts:
    tt.append(sum([pxy[i]*pxy[i+t] + pxz[i]*pxz[i+t] + pyz[i]*pyz[i+t] for i in range(0,N-t,t)])/(N-t)*t)
import pylab as pl
pl.plot(ts,tt)
pl.show()
