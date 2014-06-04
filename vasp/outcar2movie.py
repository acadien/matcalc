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
import msd

mlab.options.offscreen = True

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
mlab.figure(bgcolor=(0.8,0.8,0.8))
for i,b in enumerate(bytenums):
    if i%2==0:
        continue
    outcar.seek(b)
    outcar.readline()
    outcar.readline()
    
    atoms = [map(float,outcar.readline().split()[:3]) for a in range(nAtoms)]
    tet,rcut = orderParam.tetrahedral(atoms,basis,rcut=3.2)
    
    ax,ay,az=zip(*atoms)
    mlab.points3d(ax,ay,az,tet,colormap='jet',scale_factor=1.5,scale_mode='none',resolution=8)
    #mlab.plot3d([0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]],opacity=0)
    mlab.view(focalpoint=(basis[0][0]/1.1,basis[1][1]/1.1,basis[2][2]/1.1))
    n = str(i).zfill(5)
    mlab.savefig(movieDir+prefix+n+suffix,magnification=1.5)
    mlab.clf()

#mencoder "mf://*.png" -mf fps=25 -o output.avi -ovc lavc -lavcopts vcodec=mpeg4

#opt="vbitrate=2160000:mbd=2:keyint=132:vqblur=1.0:cmp=2:subcmp=2:dia=2:mv0:last_pred=3"
#mencoder -ovc lavc -lavcopts vcodec=msmpeg4v2:vpass=1:$opt -mf type=png:fps=48 -nosound -o /dev/null mf://\*.png
#mencoder -ovc lavc -lavcopts vcodec=msmpeg4v2:vpass=2:$opt -mf type=png:fps=48 -nosound -o output.avi mf://\*.png
