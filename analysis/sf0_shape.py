#!/usr/bin/python

#mine
import utils
import plotRemote as pr
from rdfExtend import rdfExtend
from sf import sfq0
#theirs
import sys
import pylab as pl
from math import pi,sqrt
from scipy import integrate
from datatools import superSmooth
import numpy as np
import itertools
from colors import vizSpec

utils.usage(["<rdf file, typically average>"],1,1)

rdfFile = sys.argv[1]
rdfData = open(rdfFile).readlines()[1:]

parse = lambda x: map(float,x.split())
combine = lambda x: " ".join(map(str,i))+"\n"
ndens = 288/(19.058221406**3)
rdfX,rdfY = map(np.array,zip(*map(parse,rdfData)))
rdfY = 4*pi*ndens*(rdfY-1)*rdfX

pl.plot(rdfX,rdfY)
pl.show()

dx = rdfX[1]-rdfX[0]

rmax = 100.0
T = 800.0
eps = 1.0
rm = 2.5
Niter = 9
damped = 0.04
#cr,rdfExtX,rdfExtY = rdfExtend(rdfX,rdfY,ndens,rmax=rmax,Niter=Niter,T=T,rm=rm,eps=eps,damped=damped,PY=True)

#rdfX=rdfExtX
#rdfY=rdfExtY
Lmax=19.058221406
qbins=1000
minq,maxq,dq=0,Lmax,Lmax/qbins
qs=np.array([i*dq+minq for i in range(qbins)])
qs[0]=1E-10
rdfX = np.array(rdfX)
rdfY = np.array(rdfY)

#lmn = range(100)
#qlmn = itertools.product(lmn,repeat=3)
#qs = [2*pi/Lmax*sqrt(i[0]*i[0]+i[1]*i[1]+i[2]*i[2]) for i in qlmn]
#qs.sort()
#qs = [q for q in qs if q<8.0]
#print len(qs)

def integ(x,y):
    if len(x)==1:
        return 0
    a=integrate.simps(y,x=x)
    return a

def alpha(lrs,br):
    a=list()
    #little-r, big-R
    for lr in lrs:
        if lr > 2*br:
            a.append(0.0)
        else:
            a.append((1.0 - lr/(2*br))**2 * (1.0 + lr/(4*br)))
    return np.array(a)

bigR=1000.0
sf=list()
a = alpha(rdfX,bigR)
    
for q in qs:
    qr=q*rdfX
        #pl.plot(rdfY*np.sin(qr)/q*a ,label=str(q),c=vizSpec(float(q)/max(qs)))
    sf.append(1 + integ(rdfX, rdfY*np.sin(qr)/q*a))
pl.plot(qs,sf,label=str(bigR))

pl.legend(loc=0)
pl.show()
