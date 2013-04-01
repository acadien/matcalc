#!/usr/bin/python

import sys
#mine
import plotRemoteMaya as prm

from colors import float2rgb
from scipy import array
from mayavi import mlab
from numpy import linspace,array
import matplotlib.cm as cm
#mine
import poscarIO

def usage():
    print "%s <POSCAR file1> <POSCAR file2>"%sys.argv[0]
    print "plots the vectors of atomic positions from file1 to file2."
    print "Common Example:"
    print "    poscarPlotDiff.py POSCAR CONTCAR"

if len(sys.argv) != 3:
    usage()
    exit(0)

poscar1=open(sys.argv[1],"r").readlines()
poscar2=open(sys.argv[2],"r").readlines()

[basis,atypes,atoms,head,poscar] = poscarIO.read(poscar1)
ax1,ay1,az1=map(array,zip(*atoms))

[basis,atypes,atoms,head,poscar] = poscarIO.read(poscar2)
ax2,ay2,az2=map(array,zip(*atoms))

dx = ax1-ax2
dy = ay1-ay2
dz = az1-az2

fig=mlab.figure(bgcolor=(1.0,1.0,1.0))

d=((dx**2+dy**2+dz**2)/3)**0.5

mlab.quiver3d(ax1,ay1,az1,dx,dy,dz,scalars=d,scale_factor=1.0,scale_mode='none',line_width=4.0)
cb = mlab.colorbar(title="Distance in Angstroms",)
cb._title_text_property.color = (0.0,0.0,0.0)
cb._label_text_property.color = (0.0,0.0,0.0)
mlab.show()
