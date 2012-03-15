#!/usr/bin/python

import sys
import pylab as pl
from scipy import array
from mpl_toolkits.mplot3d import Axes3D
from math import *
#mine
from poscarIO import readposcar
from voronoiNeighbors import voronoiNeighbors
from struct_tools import *
from colors import float2rgb
from paircor import paircor_ang

def usage():
    print "%s <POSCAR file> <bond length,err>"%sys.argv[0]

if len(sys.argv) < 2:
    usage()
    exit(0)

bl=-1.
bw=-1.
if len(sys.argv) >= 3:
    bl,bw=map(float,sys.argv[2].split(","))

poscar=open(sys.argv[1],"r").readlines()

[v1,v2,v3,atypes,ax,ay,az,head,poscar] = readposcar(poscar)
atoms=array(zip(ax,ay,az))
basis=array([v1,v2,v3])
lengths=array([v1[0],v2[1],v3[2]])

j=0
types=list()
for i in atypes:
    types+=[j]*i
    j+=1

neighbs=voronoiNeighbors(atoms=atoms,atypes=atypes,basis=basis,style='full')

fig=pl.figure()
aa = fig.add_subplot(111,projection='3d')
aa.scatter(*zip(*atoms),c='blue')
specbonds=[list() for i in range(len(atoms))]
for i in range(len(atoms)):
    ai=atoms[i]
    for j in neighbs[i]:
        aj=array(list(atoms[j]))
        if dist(ai,aj)!=dist_periodic(ai,aj,lengths):
            for c in range(3):
                 d = aj[c]-ai[c]
                 if d>lengths[c]/2.0: aj[c] -= lengths[c]
                 elif d<-lengths[c]/2.0: aj[c] += lengths[c]
        d=dist(ai,aj)
        c=float2rgb(d,3.,6.)
        ls='.-'
        if fabs(d-bl)<bw:
            c='black'
            ls='-'
            specbonds[i].append(j)
        aa.plot(*zip(ai,aj),c=c,ls=ls)
pl.title(sys.argv[1])


#Analyze specbond lengths
xb,yb=paircor_ang(atoms,specbonds,basis)
pl.figure()
pl.plot(xb,yb)

pl.show()    
    
