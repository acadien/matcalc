#!/usr/bin/python

import sys
#mine
import plotRemoteMaya as prm

import poscarIO
from voronoiNeighbors import voronoiNeighbors
from struct_tools import *
from colors import float2rgb

from scipy import array
from math import *
from mayavi import mlab


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

[basis,atypes,atoms,head,poscar] = poscarIO.read(poscar)
ax,ay,az=zip(*atoms)
v1,v2,v3=basis
j=0
types=list()
for i in atypes:
    types+=[j+1]*i
    j+=1

fig=mlab.figure()

mlab.points3d(ax,ay,az,types,scale_factor=1.0,vmin=-1.0,vmax=4.0)
z=[0,0,0]
mlab.plot3d([0,v1[0]],[0,v1[1]],[0,v1[2]],color=(1,1,1),line_width=0.1)
mlab.plot3d([0,v2[0]],[0,v2[1]],[0,v2[2]],color=(1,1,1),line_width=0.1)
mlab.plot3d([0,v3[0]],[0,v3[1]],[0,v3[2]],color=(1,1,1),line_width=0.1)

prm.prmshow(fname="POSCAR.png")
