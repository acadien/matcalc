#!/usr/bin/python

import sys
from math import *
from scipy import array,zeros
import pylab as pl
import subprocess
import os
from mayavi import mlab
#mine
import outcarIO
import orderParam

#make a directory for the movie
#for configs in outcar
#  read in a config from outcar
#  calculate tetrahedral order param or CN
#  generate the image and save as a file

def usage():
    print "Usage: %s <Outcar> <Reference Config #>"%sys.argv[0].split("/")[-1]

if len(sys.argv)<2:
    usage()
    exit(0)
    
outcarFile = sys.argv[1]
nAtoms = outcarIO.nIons(outcarFile)
basis = array(map(array,outcarIO.basis(outcarFile)))

movieDir = "./outcarMovie/"
prefix = "outcarTet"
suffix = ".png"
if not os.path.exists(movieDir): os.makedirs(movieDir)

#Find the starting locations of atomic data in outcarfile
grepResults = subprocess.check_output("grep -b POSITION %s"%outcarFile,shell=True).split("\n")
bytenums=[int(i.split(":")[0]) for i in grepResults if len(i)>2]

outcar= open(outcarFile,"r")

for i,b in enumerate(bytenums):
    outcar.seek(b)
    outcar.readline()
    outcar.readline()
    
    atoms = [map(float,outcar.readline().split()[:3]) for a in range(nAtoms)]
    tet,rcut = orderParam.tetrahedral(atoms,basis,rcut=3.2)
    
    ax,ay,az=zip(*atoms)
    mlab.points3d(ax,ay,az,tet,colormap='jet',scale_factor=1.9,scale_mode='none',resolution=8)
    
    n = str(i).zfill(5)
    mlab.savefig(movieDir+prefix+n+suffix,magnification=1.5)
    mlab.clf()

#mencoder "mf://*.png" -mf fps=25 -o output.avi -ovc lavc -lavcopts vcodec=mpeg4
