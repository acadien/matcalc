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

#vhalfNeighbors=voronoiNeighbors(atoms=atoms,basis=basis,atypes=atypes,style='half')
halfNeighbors=neighbors(atoms,array([[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]),max([l+2*w for l,w in zip(lengths,widths)]),style='half')

fullNeighbors=half2full(halfNeighbors)
coordNumbers=map(len,fullNeighbors)
print coordNumbers
for l,w in zip(lengths,widths):
    atomPairs=atomsAtLength(atoms,halfNeighbors,l,w,periodic=True,bounds=bounds)

    CNhist=zeros(100)
    cns=list()
    for [a1,a2] in atomPairs:
        c1=coordNumbers[a1]
        c2=coordNumbers[a2]
        cns.append(c1)
        cns.append(c2)
        CNhist[c1]+=1
        CNhist[c2]+=2

        #pl.plot(CNhist,label=str(l)+"$\pm$"+str(w))
    mn=max(min(cns)-5,0)
    mx=min(max(cns)+5,len(CNhist))
                         
    pl.plot(range(mn,mx+1),CNhist[mn:mx+1],label=str(l)+"$\pm$"+str(w))
pl.xlabel("Coordination Number")
#pl.xlim([i
pl.ylabel("Count")
pl.legend(title="Atoms with Bond Length")
pl.show()
