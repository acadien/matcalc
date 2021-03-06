#!/usr/bin/python

import sys
from math import *
import pylab as pl
from scipy import *
#mine
import poscarIO
from geometry import atomsAtLength
from neighbors import voronoiNeighbors,half2full
from struct_tools import *

if len(sys.argv)<2:
    print "Usage:"
    print sys.argv[0]+" <poscar> <opt: max bond length>"
    exit(0)

poscar = open(sys.argv[1],"r").readlines()
mbl=-1.
if len(sys.argv)==3:
    mbl=float(sys.argv[2])

[basis,atypes,atoms,head,poscar] = poscarIO.read(poscar)
bounds=array([basis[0][0],basis[1][1],basis[2][2]])

#Generate neighbor lists, half and full
vhalfNeighbors=voronoiNeighbors(atoms,basis,atypes,style='half')
vfullNeighbors=half2full(vhalfNeighbors)

#Check lengths of bonds
if mbl>0:
    ls=list()
    for ai,neighbs in enumerate(vfullNeighbors):
        ls.append([dist_periodic(atoms[ai],atoms[aj],bounds) for aj in neighbs])
        
    vfullNeighbors= [[vfullNeighbors[i][j] for j,dist in enumerate(dists) if dist<=mbl] for i,dists in enumerate(ls)]
    for i,j in zip(vfullNeighbors,ls):
        print len(i),len(j)

coordNumbers=map(len,vfullNeighbors)

CNhist=zeros(200)
for c in coordNumbers:
    CNhist[c]+=1
mn=max(min(coordNumbers)-5,0)
mx=min(max(coordNumbers)+5,len(CNhist))
                         
pl.plot(range(mn,mx+1),CNhist[mn:mx+1])

pl.xlabel("Coordination Number")
pl.ylabel("Count")
pl.show()
