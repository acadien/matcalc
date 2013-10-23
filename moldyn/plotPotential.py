#!/usr/bin/python

#plot multiple potentials on top of eachother!
#start out with just ADP

import plotRemote as pr
import matplotlib
from math import *

import sys
import pylab as pl
from numpy import array 

def usage():
    print "%s <adp potential files>"%sys.argv[0].split("/")[-1]

def loadNparseADP(filename):
    lmpdata=open(filename,"r").readlines()

    #load up adp data
    nrho,drho,nr,dr,cut=lmpdata[4].split()
    nrho=int(nrho)
    drho=float(drho)
    nr=int(nr)
    dr=float(dr)
    cut=float(cut)
    #F(rho),phi,rho,u,w
    lmpdata=lmpdata[6:]
    lydata=[list() for i in range(5)]
    lxdata=[list() for i in range(5)]

    lxdata[0]=map(lambda x:x*drho,range(nrho)) #F
    lydata[0]=map(float,lmpdata[:nrho])
    lxdata[1]=map(lambda x:x*dr,range(nr)) #rho
    lydata[1]=map(float,lmpdata[nrho:nrho+nr])
    lxdata[2]=map(lambda x:x*dr,range(nr)) #phi
    lydata[2]=map(float,lmpdata[nrho+nr:nrho+2*nr])
    lydata[2]=[lydata[2][i]/(lxdata[2][i]+1E-20) for i in range(nr)]
    lxdata[3]=map(lambda x:x*dr,range(nr)) #u
    lydata[3]=map(float,lmpdata[nrho+2*nr:nrho+3*nr])
    lxdata[4]=map(lambda x:x*dr,range(nr)) #w
    lydata[4]=map(float,lmpdata[nrho+3*nr:nrho+4*nr])

    varnames = ["F","rho","phi","u","w"]
    return lxdata,lydata,varnames

lxs,lys=list(),list()
for filename in sys.argv[1:]:
    lx,ly,varnames = loadNparseADP(filename)
    lxs.append(lx)
    lys.append(ly)

#print "Variable : minval - maxval"
pl.figure(figsize=(9,6))
for i,var in enumerate(varnames):

    pl.subplot(321+i)

    for j in range(len(lxs)):
        pl.plot(lxs[j][i],lys[j][i],label="/\n".join(sys.argv[1+j].split("/")[-2:]))
    dummy,xmax = pl.xlim()
    ymin,ymax = pl.ylim()
    xmax*=0.85
    if ymax>0:
        ymax*=0.85
    else:
        ymax/=0.65
    pl.annotate(var,[xmax,ymax],[xmax,ymax])

pl.legend(bbox_to_anchor=[2.0,1])
    
#pl.legend(loc=0)
pr.prshow("SampledPotential.eps")
#pl.show()
