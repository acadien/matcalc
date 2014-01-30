#!/usr/bin/python

#ADD capability to plot multiple potentials on top of eachother!

import plotRemote as pr
import matplotlib
from math import *

import sys
import pylab as pl
from numpy import array 

def usage():
    print "%s <potential sampled file> <lammps potential>"%sys.argv[0].split("/")[-1]

if len(sys.argv)<3:
    usage()
    exit(0)

potdata=open(sys.argv[1],"r").readlines()
lmpdata=open(sys.argv[2],"r").readlines()

nvars=sum([1 for i in potdata if len(i)<=1])/2
if nvars==5:
    varnames=["F","rho","phi","u","w"]

datax=[[]]
datay=[[]]
for i in potdata:
    if len(i)<=1:
        if len(datax[-1])>0:
            datax.append(list())
            datay.append(list())
        continue
    x,y,dumb=i.split()
    datax[-1].append(float(x))
    datay[-1].append(float(y))
#swap rho and F
tempx=datax[0]
tempy=datay[0]
datax[0]=datax[2]
datax[2]=tempx
datay[0]=datay[2]
datay[2]=tempy


#load up lammps data
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
print lxdata[2]
#print "Variable : minval - maxval"
pl.figure(figsize=(9,6))
for i,var in enumerate(varnames):

    pl.subplot(321+i)

    pl.scatter(datax[i],datay[i])
    xlims=pl.xlim()
    ylims=pl.ylim()

    pl.plot(lxdata[i],lydata[i],label=var)
    pl.xlim(xlims)
    pl.ylim(ylims)
#    prit "%3.3s : %4.4e - %4.4e"%(var,min(datay[i]),max(datay[i]))
    pl.legend(loc=0)
    

#pl.legend()
    #pl.savefig("/home/acadien/Dropbox/stuff%s.png"%var)
pl.legend(loc=0)
pr.prshow("SampledPotential.eps")
#pl.show()
