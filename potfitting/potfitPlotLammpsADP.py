#!/usr/bin/python

#ADD capability to plot multiple potentials on top of eachother!

import plotRemote as pr
import matplotlib
from math import *

import sys
import pylab as pl
from numpy import *

def nPotentials(x):
    nPhi = x*(x+1)/2.
    nRho = x
    nF = x
    nU = x*(x+1)/2.
    nW = x*(x+1)/2.
    return map(int, [nPhi,nRho,nF,nU,nW])

def parsePfknots(pfend):
    pfData=open(pfend,"r").readlines()

    #Parse the header
    dummy, nt, n = pfData[0].split()
    if nt!='3':
        print "This parser can only handle type 3 potfit formats for ADP potentials"
        exit(0)\

    while True:
        line = pfData.pop(0)
        if line[:2] == "#C":
            elems=line[2:].strip().split()
        if line[:2] == "#G":
            n = len(line.split())-1
            nElement = (-7+sqrt(49+24*float(n)))/6
        if line[:2] == "#E":
            break
    print "nElement types :",nElement
    nPots = nPotentials(nElement)
    xPhis,yPhis,yPhiDerivs = list(),list(),list()
    xRhos,yRhos,yRhoDerivs = list(),list(),list()
    xFs,yFs,yFDerivs = list(),list(),list()
    xUs,yUs,yUDerivs = list(),list(),list()
    xWs,yWs,yWDerivs = list(),list(),list()

    #Drop the spaces and comments
    pfData = [line for line in pfData if len(line.strip())>0 and line[0]!="#"]

    #X-Data
    genXs=lambda x:arange(x[0],x[1]+0.1,(x[1]-x[0])/(x[2]-1))
    for i in range(nPots[0]):
        l=map(float,pfData.pop(0).split())
        xPhis.append(genXs(l))

    for i in range(nPots[1]):
        l=map(float,pfData.pop(0).split())
        xRhos.append(genXs(l))

    for i in range(nPots[2]):
        l=map(float,pfData.pop(0).split())
        xFs.append(genXs(l))

    for i in range(nPots[3]):
        l=map(float,pfData.pop(0).split())
        xUs.append(genXs(l))

    for i in range(nPots[4]):
        l=map(float,pfData.pop(0).split())
        xWs.append(genXs(l))

    #Y-data
    for xPhi in xPhis:
        yPhiDerivs.append(map(float,pfData[0].split()))
        yPhis.append(array(map(float,pfData[1:size(xPhi)+1])))
        pfData=pfData[1+size(xPhi):]

    for xRho in xRhos:
        yRhoDerivs.append(map(float,pfData[0].split()))
        yRhos.append(array(map(float,pfData[1:size(xRho)+1])))
        pfData=pfData[1+size(xRho):]

    for xF in xFs:
        yFDerivs.append(map(float,pfData[0].split()))
        yFs.append(array(map(float,pfData[1:size(xF)+1])))
        pfData=pfData[1+size(xF):]

    for xU in xUs:
        yUDerivs.append(map(float,pfData[0].split()))
        yUs.append(array(map(float,pfData[1:size(xU)+1])))
        pfData=pfData[1+size(xU):]

    for xW in xWs:
        yWDerivs.append(map(float,pfData[0].split()))
        yWs.append(array(map(float,pfData[1:size(xW)+1])))
        pfData=pfData[1+size(xW):]

    #Insert new points at both ends to account for derivatives
    step=0.1
    def appendDerivs(x,y,primes,step):
        x=insert(x,0,[x[0]-step])
        y=insert(y,0,[y[0]-primes[0]*step])

        x=insert(x,-1,[x[-1]-step])
        y=insert(y,-1,[y[-1]-primes[1]*step])
        return x,y

    xPhis,yPhis = zip(*[appendDerivs(xPhi,yPhi,yPhiDeriv,step) \
                  for xPhi,yPhi,yPhiDeriv in zip(xPhis,yPhis,yPhiDerivs)])
        
    xRhos,yRhos = zip(*[appendDerivs(xRho,yRho,yRhoDeriv,step) \
                  for xRho,yRho,yRhoDeriv in zip(xRhos,yRhos,yRhoDerivs)])

    xFs,yFs = zip(*[appendDerivs(xF,yF,yFDeriv,step) \
              for xF,yF,yFDeriv in zip(xFs,yFs,yFDerivs)])

    xUs,yUs = zip(*[appendDerivs(xU,yU,yUDeriv,step) \
              for xU,yU,yUDeriv in zip(xUs,yUs,yUDerivs)])

    xWs,yWs = zip(*[appendDerivs(xW,yW,yWDeriv,step) \
              for xW,yW,yWDeriv in zip(xWs,yWs,yWDerivs)])

    return [xFs,xRhos,xPhis,xUs,xWs],[yFs,yRhos,yPhis,yUs,yWs] #return in lammps order!

def parseLammpsADP(fname):
    lmpdata=open(fname,"r").readlines()

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
    return lxdata,lydata

def usage():
    print "%s <lammps potential> <optional pf_end>"%sys.argv[0].split("/")[-1]

if len(sys.argv)<2:
    usage()
    exit(0)

lammpsx,lammpsy =parseLammpsADP(sys.argv[1])

pf=False
if len(sys.argv)==3:
    pf=True
    pfknotsx,pfknotsy = parsePfknots(sys.argv[2])

varnames=["F","rho","phi","u","w"]

#print "Variable : minval - maxval"
pl.figure(figsize=(9,6))
for i,var in enumerate(varnames):

    pl.subplot(321+i)

    if pf:
        pl.scatter(pfknotsx[i],pfknotsy[i])
        xlims=pl.xlim()
        ylims=pl.ylim()

    pl.plot(lammpsx[i],lammpsy[i],label=var)

    if pf:
        pl.xlim(xlims)
        pl.ylim(ylims)
        if i==2:
            pl.ylim([0,2])

    pl.legend(loc=0)
    

#pl.legend()
    #pl.savefig("/home/acadien/Dropbox/stuff%s.png"%var)
pl.legend(loc=0)
pr.prshow("SampledPotential.eps")
#pl.show()
