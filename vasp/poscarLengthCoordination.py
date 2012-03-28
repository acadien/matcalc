#!/usr/bin/python

import sys
from math import *
import pylab as pl
from scipy import *
#mine
from poscarIO import readposcar
from geometry import atomsAtLength
from voronoiNeighbors import voronoiNeighbors
from struct_tools import *

if len(sys.argv)<3:
    print "Usage:"
    print sys.argv[0]+" <poscar> <length0>,<width0> <length1>,<width1> ..."
    exit(0)

poscar = open(sys.argv[1],"r").readlines()
lengths,widths = zip(*[map(float,i.split(",")) for i in sys.argv[2:]])

[v1,v2,v3,atypes,ax,ay,az,head,poscar] = readposcar(poscar)
atoms=array(zip(ax,ay,az))
basis=[v1,v2,v3]
bounds=array([v1[0],v2[1],v3[2]])

vhalfNeighbors=voronoiNeighbors(atoms=atoms,basis=basis,atypes=atypes,style='half')
halfNeighbors=neighbors(atoms,array([[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]),max([l+w for l,w in zip(lengths,widths)]),style='half')

vfullNeighbors=half2full(vhalfNeighbors)
fullNeighbors=half2full(halfNeighbors)

print set(flatten(fullNeighbors))

vCoordNumbers=map(len,vfullNeighbors)
coordNumbers=map(len,fullNeighbors)

for l,w in zip(lengths,widths):
    atomPairs=atomsAtLength(atoms,halfNeighbors,l,w,periodic=True,bounds=bounds) 
    vAtomPairs=atomsAtLength(atoms,vhalfNeighbors,l,w,periodic=True,bounds=bounds)

    #For the bond length in question, how many bonds of length in question?
    CNblqhist=zeros(100)
    cnsblq=list()
    for [a1,a2] in atomPairs:
        c1=coordNumbers[a1]
        c2=coordNumbers[a2]
        cnsblq.append(c1)
        cnsblq.append(c2)
        CNblqhist[c1]+=1
        CNblqhist[c2]+=1
    mnblq=max(min(cnsblq)-3,0)
    mxblq=min(max(cnsblq)+3,len(CNblqhist))                   
    pl.plot(range(mnblq,mxblq+1),CNblqhist[mnblq:mxblq+1],label=str(l)+"$\pm$"+str(w))


    #For atoms with a bond of length in question, how many total bonds?
    CNtothist=zeros(100)
    cnstot=list()
    for [a1,a2] in vAtomPairs:
        c1=vCoordNumbers[a1]
        c2=vCoordNumbers[a2]
        cnstot.append(c1)
        cnstot.append(c2)
        CNtothist[c1]+=1
        CNtothist[c2]+=1
    mntot=max(min(cnstot)-3,0)
    mxtot=min(max(cnstot)+3,len(CNtothist))                   
    pl.plot(range(mntot,mxtot+1),CNtothist[mntot:mxtot+1],label=str(l)+"$\pm$"+str(w))


pl.legend(["Atoms with Bond Length, how many bonds of length in question.","Atoms with Bond Length, how many total bonds."])
pl.xlabel("Coordination Number")
pl.ylabel("Count")
pl.show()


