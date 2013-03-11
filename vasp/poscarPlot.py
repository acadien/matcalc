#!/usr/bin/python

import sys
#mine
import plotRemoteMaya as prm

import poscarIO
from colors import float2rgb
from scipy import array
from mayavi import mlab
import matplotlib.cm as cm
#mine
#from voronoiNeighbors import voronoiNeighbors
#from struct_tools import *
from orderParam import coordinationNumbers

orderParams={"cn":coordinationNumbers}

def usage():
    print "%s <POSCAR file> <order parameter>"%sys.argv[0]
    print "Order Parameter can be one of:"
    print "   cn : Coordination Number"
    print ""

if len(sys.argv) < 2:
    usage()
    exit(0)

opFlag = False
if len(sys.argv)==3:
    op = sys.argv[2]
    opFlag = True

poscar=open(sys.argv[1],"r").readlines()

[basis,atypes,atoms,head,poscar] = poscarIO.read(poscar)
ax,ay,az=zip(*atoms)
v1,v2,v3=basis
j=0
types=list()
for i in atypes:
    types+=[j+1]*i
    j+=1

if opFlag:
    ops = orderParams[op](array(atoms),array(basis))
    if max(ops)==min(ops):
        opsn=[1]*len(ops)
    else:
        opsn = [(i-float(min(ops)))/float(max(ops)-min(ops)) for i in ops]
fig=mlab.figure()

if opFlag:
    reorg={}
    for o in set(opsn):
        reorg[o]=[v for j,v in enumerate(atoms) if opsn[j]==o]  
    for i in reorg.keys():
        aa=reorg[i]
        ax,ay,az=zip(*aa)
        mlab.points3d(ax,ay,az,color=cm.jet(i)[:3],scale_factor=1.0,vmin=0.0,vmax=1.0)
    print sorted(set(ops))
else:
    mlab.points3d(ax,ay,az,ops,scale_factor=1.0,vmin=-1.0,vmax=4.0)
z=[0,0,0]
mlab.plot3d([0,v1[0]],[0,v1[1]],[0,v1[2]],color=(1,1,1),line_width=0.1)
mlab.plot3d([0,v2[0]],[0,v2[1]],[0,v2[2]],color=(1,1,1),line_width=0.1)
mlab.plot3d([0,v3[0]],[0,v3[1]],[0,v3[2]],color=(1,1,1),line_width=0.1)

prm.prmshow(fname="POSCAR.png")
