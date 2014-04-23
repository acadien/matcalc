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
from meanSquareDist import meanSquareDist
import poscarIO

orderParams={"CN":coordinationNumber, \
             "BO":bondOrientation, \
             "TET":tetrahedral,\
             "MSD":meanSquareDist    }

def usage():
    print "%s <order parameter> <POSCAR file>"%sys.argv[0]
    print "Order Parameter can be one of:"
    print "     CN : Coordination Number"
    print "    BO# : Bond Orientation (Q) with l=#"
    print "    TET : Tetrahedral order parameter (Sg)"
    print "Flags can be used anywhere in args:"
    print "-rectify: applies PBC to atoms on a cubic cell"
    print "   -rcut: cutoff distance when building neighbor list, follow by float value"
    print ""

if len(sys.argv) < 2:
    usage()
    exit(0)

rectifyFlag=False
rcut=None
opFlag = False
lval=0
for i,v in enumerate(sys.argv):
    if v in ["-rectify","-Rectify"]:
        rectifyFlag=True
        sys.argv.pop(i)
    if v == "-rcut":
        rcut=float(sys.argv[i+1])
        sys.argv.pop(i) #pop the flag
        sys.argv.pop(i) #pop the float

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

if rectifyFlag:
    atoms=[[i[0]%basis[0][0],i[1]%basis[1][1],i[2]%basis[2][2]]for i in atoms]

ax,ay,az=zip(*atoms)
v1,v2,v3=basis
j=0
types=list()
for i in atypes:
    types+=[j+1]*i
    j+=1

fig=mlab.figure()#bgcolor=(1.0,1.0,1.0))

#Get the order parameter and convert to integer format (opsn) for
#coloring of atoms
if opFlag:
    ops,rcut = orderParams[op](array(atoms),array(basis),l=lval,rcut=rcut)
    print len(ops)
    n=None
    mnop = min(ops)
    mxop = max(ops)
    #if op=="TET":
        #mnop=0.0
        #mxop=1.0
        #n=11

    #Only plot with coloring if there is some variance in the OP
    if mxop - mnop > 1E-10:
        #spectral
        mp3d = mlab.points3d(ax,ay,az,ops,colormap='jet',scale_factor=1.9,scale_mode='none',resolution=14)
        if n==None:
            n=min(len(set(ops)),10)
    else:
        mp3d = mlab.points3d(ax,ay,az,[mnop]*len(az),colormap='jet',scale_factor=1.9,scale_mode='none',resolution=14)
        n=2

    #Color bar formatting
    cb = mlab.colorbar(title=sys.argv[1], orientation='vertical', nb_labels=n,nb_colors=n)
    cb.use_default_range = False
    cb.data_range = (mnop,mxop)
    
else:
    mp3d = mlab.points3d(ax,ay,az,types,scale_factor=2.0,scale_mode='none',color=(0.25,0.33,1.00),resolution=14)

z=[0,0,0]
#extent = (0,v1[0],0,v2[0],
#mlab.outline(mp3d, color=(.7, .7, .7), extent=cat1_extent)
mlab.plot3d([0,v1[0],v1[0]+v2[0],v2[0],0,v3[0]],[0,v1[1],v1[1]+v2[1],v2[1],0,v3[1]],[0,v1[2],v1[2]+v2[2],v2[2],0,v3[2]],color=(0,0,0),line_width=0.5)
mlab.plot3d([v3[0],v3[0]+v1[0],v3[0]+v2[0]+v1[0],v3[0]+v2[0],v3[0]],[v3[1],v3[1]+v1[1],v3[1]+v2[1]+v1[1],v3[1]+v2[1],v3[1]],[v3[2],v3[2]+v1[2],v3[2]+v2[2]+v1[2],v3[2]+v2[2],v3[2]],color=(0,0,0),line_width=0.5)
mlab.plot3d([v1[0],v1[0]+v3[0]],[v1[1],v1[1]+v3[1]],[v1[2],v1[2]+v3[2]],color=(0,0,0),line_width=0.5)
mlab.plot3d([v2[0],v2[0]+v3[0]],[v2[1],v2[1]+v3[1]],[v2[2],v2[2]+v3[2]],color=(0,0,0),line_width=0.5)
mlab.plot3d([v1[0]+v2[0],v1[0]+v2[0]+v3[0]],[v1[1]+v2[1],v1[1]+v2[1]+v3[1]],[v1[2]+v2[2],v1[2]+v2[2]+v3[2]],color=(0,0,0),line_width=0.5)
#mlab.plot3d([0,v2[0]],[0,v2[1]],[0,v2[2]],color=(0,0,0),line_width=0.5)
#mlab.plot3d([0,v3[0]],[0,v3[1]],[0,v3[2]],color=(0,0,0),line_width=0.5)
prm.prmshow(fname="POSCAR.png")
