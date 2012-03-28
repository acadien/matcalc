#!/usr/bin/python

#The CHGCARs analyzed by this script need to be part of rectangular simulations
#parallepiped simulations can be changed to cubic/rectangular using the script
#poscarRectify.py
#This greatly simplifies the application of periodic boundary conditions

import sys
import operator
from scipy import *
import pylab as pl
#mine
from struct_tools import dist_periodic
from chgcarIO import readchgcar
from voronoiNeighbors import voronoiNeighbors
from fieldPointAnalysis import fieldNeighbors

def usage():
    print "Usage:"
    print "%s <CHGCAR> <minL0,maxL0> <minL1,maxL1>"%(sys.argv[0])
    print "Note, CHGCAR must be part of a cubic/rectangular simulation"

if len(sys.argv)<3:
    usage()
    exit(0)

BLs= [map(float,i.split(",")) for i in sys.argv[2:]]
nBLs=len(BLs)

#Parse CHGCAR
chgcar=open(sys.argv[1],"r").readlines()
(v1,v2,v3,atypes,axs,ays,azs,header),gridSize,chg = readchgcar(chgcar)
basis=asarray([v1,v2,v3])
lengths=array([v1[0],v2[1],v3[2]])
atoms=array(zip(axs,ays,azs))

#Grid properties
nGridPoints=reduce(operator.mul,gridSize)

#Neighbors
halfNeighbors=voronoiNeighbors(atoms=atoms,basis=basis,atypes=atypes,style='half')

a=chg.shape
print "Average CHG value:",sum([sum([sum(line)/a[2] for line in plane])/a[1] for plane in chg])/a[0]

#Evaluate the CHG between each nieghbor pair
avgchgline,xchglines,ychglines,halfNeighbors=fieldNeighbors(atoms,atypes,basis,chg,gridSize,halfNeighbors)

#Cutoff neighbors that fall below the thresh-hold
Ninterp=len(avgchgline)
avgIBL=zeros([nBLs,Ninterp]) #avg line inside the bond length
avgOBL=zeros([nBLs,Ninterp]) #avg line outside the bond length
nibl=zeros(nBLs)
nobl=zeros(nBLs)

cnt=0
pl.figure()
for a,jNeighbors in enumerate(halfNeighbors):
    atoma=atoms[a]
    for b,jNeighb in enumerate(jNeighbors):
        atomb=atoms[jNeighb]
        d=dist_periodic(atoma,atomb,lengths)
        vals=ychglines[a][b]
        xx=xchglines[a][b]
        for j,bl in enumerate(BLs):
            [blmin,blmax]=bl
            if d>=blmin and d<blmax:
                pl.plot(xx,vals)
                nibl[j]+=1
                avgIBL[j]+=array(vals)                   
            else:
                cnt+=1
                nobl[j]+=1
                avgOBL[j]+=array(vals)                   

for i in range(nBLs):
    if nibl[i]==0:
        nibl[i]=1
    if nobl[i]==0:
        nobl[i]=1
avgOBL=[avgOBL[i]/nobl[i] for i in range(nBLs)]
avgIBL=[avgIBL[i]/nibl[i] for i in range(nBLs)]
fige=pl.figure()
#sml=fige.add_subplot(212)
big=fige.add_subplot(111)

for i,(avgo,avgi) in enumerate(zip(list(avgOBL),list(avgIBL))):
    big.plot(avgi,label="-".join(map(str,BLs[i]))+" #B="+str(int(nibl[i])))
    #sml.plot(avgo)

#pl.legend(["Good Neighbs","Bad Neighbs"])
#sml.set_title("Average CHG Density for Inside Bonds")
#big.legend(title="Bond Len(nm)",bbox_to_anchor=(1.10,1,))
big.legend(title="Bond Len (nm)",loc=0)
#big.set_title("Average CHG Density for Outside Bonds")
big.set_title(sys.argv[1])

pl.show()

