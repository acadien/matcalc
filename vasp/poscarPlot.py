#!/usr/bin/python

import sys
#mine
import plotRemoteMaya as prm

from colors import float2rgb
from scipy import array
from mayavi import mlab
from numpy import linspace
import matplotlib.cm as cm
#mine
from orderParam import coordinationNumber,bondOrientation,tetrahedral
import poscarIO

orderParams={"CN":coordinationNumber, \
             "BO":bondOrientation, \
             "TET":tetrahedral}

def usage():
    print "%s <order parameter> <POSCAR file>"%sys.argv[0]
    print "Order Parameter can be one of:"
    print "   CN : Coordination Number"
    print "   BO# : Bond Orientation (Q) with l=#"
    print "   TET : Tetrahedral order parameter (Sg)"
    print ""

if len(sys.argv) < 2:
    usage()
    exit(0)

opFlag = False
lval=0
if len(sys.argv)==2:
    poscar=open(sys.argv[1],"r").readlines()
elif len(sys.argv)==3:
    op = sys.argv[1]
    if op[:2]=="BO":
        lval=int(op[-1])
        op="BO"
    opFlag = True
    poscar=open(sys.argv[2],"r").readlines()

[basis,atypes,atoms,head,poscar] = poscarIO.read(poscar)
ax,ay,az=zip(*atoms)
v1,v2,v3=basis
j=0
types=list()
for i in atypes:
    types+=[j+1]*i
    j+=1

fig=mlab.figure(bgcolor=(1.0,1.0,1.0))

#Get the order parameter and convert to integer format (opsn) for
#coloring of atoms
if opFlag:
    ops = orderParams[op](array(atoms),array(basis),lval)
    mnop = min(ops)
    mxop = max(ops)

    #Only plot with coloring if there is some variance in the OP
    if mxop - mnop > 1E-10:
        mlab.points3d(ax,ay,az,ops,colormap='jet',scale_factor=1.0,scale_mode='none')
        n=min(len(set(ops)),10)
    else:
        mlab.points3d(ax,ay,az,[mnop]*len(az),scale_factor=1.0,scale_mode='none')
        n=2

    #Color bar formatting
    cb = mlab.colorbar(title=sys.argv[1], orientation='vertical', nb_labels=n,nb_colors=n)
    cb.use_default_range = False
    cb.data_range = (min(ops),max(ops))
    
else:
    mlab.points3d(ax,ay,az,types,scale_factor=1.0,scale_mode='none')

z=[0,0,0]
mlab.plot3d([0,v1[0]],[0,v1[1]],[0,v1[2]],color=(0,0,0),line_width=0.5)
mlab.plot3d([0,v2[0]],[0,v2[1]],[0,v2[2]],color=(0,0,0),line_width=0.5)
mlab.plot3d([0,v3[0]],[0,v3[1]],[0,v3[2]],color=(0,0,0),line_width=0.5)

prm.prmshow(fname="POSCAR.png")
