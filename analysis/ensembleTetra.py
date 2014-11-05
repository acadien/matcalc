#!/usr/bin/python

import sys
from math import *
import subprocess
import os
#mine
import outcarIO
import lammpsIO
import orderParam

from scipy import array,zeros
import pylab as pl

def usage():
    print "Usage: %s <Outcar> "%sys.argv[0].split("/")[-1]

if len(sys.argv)<2:
    usage()
    exit(0)

rcut = 3.2

filename = sys.argv[1]
outcarFlag=False
lammpsFlag=False
if "OUTCAR" in filename:
    outcarFlag=True
    lammpsFlag=False
else:
    outcarFlag=False
    lammpsFlag=True


if outcarFlag:
    nAtoms = outcarIO.nIons(filename)
    basis = array(map(array,outcarIO.basis(filename)))

    #Find the starting locations of atomic data in outcarfile
    grepResults = subprocess.check_output("grep -b POSITION %s"%filename,shell=True).split("\n")
    bytenums=[int(i.split(":")[0]) for i in grepResults if len(i)>2]

    outcar= open(filename,"r")
    tetraOut =open(filename+".tetra","w")

    tetraOut.write("AverageTetra PerAtomTetra\n")
    for i,b in enumerate(bytenums):
        outcar.seek(b)
        outcar.readline()
        outcar.readline()

        atoms = [map(float,outcar.readline().split()[:3]) for a in range(nAtoms)]
        tet,rcut = orderParam.tetrahedral(atoms,basis,rcut=rcut)
        tAvg = sum(tet)/len(tet)
        tets = " ".join(map(str,tet))
        tetraOut.write(str(tAvg)+" "+tets+"\n")


if lammpsFlag:
    nAtoms = lammpsIO.nAtoms(filename)
    basisByteNums = lammpsIO.basisBytes(filename)
    atomsByteNums = lammpsIO.atomsBytes(filename)
    
    tetraOut = open(filename+".tetra","w")

    tetraOut.write("AverageTetra PerAtomTetra\n")
    for i,(bByte,aByte) in enumerate(zip(basisByteNums,atomsByteNums)):
        basis = lammpsIO.parseBasis(filename,bByte)

        atoms,dummy = lammpsIO.parseAtoms(filename,aByte,nAtoms,basis)

        tet,rcut = orderParam.tetrahedral(atoms,basis,rcut=rcut)
        tAvg = sum(tet)/len(tet)
        tets = " ".join(map(str,tet))
        tetraOut.write(str(tAvg)+" "+tets+"\n")

#mencoder "mf://*.png" -mf fps=25 -o output.avi -ovc lavc -lavcopts vcodec=mpeg4

#opt="vbitrate=2160000:mbd=2:keyint=132:vqblur=1.0:cmp=2:subcmp=2:dia=2:mv0:last_pred=3"
#mencoder -ovc lavc -lavcopts vcodec=msmpeg4v2:vpass=1:$opt -mf type=png:fps=48 -nosound -o /dev/null mf://\*.png
#mencoder -ovc lavc -lavcopts vcodec=msmpeg4v2:vpass=2:$opt -mf type=png:fps=48 -nosound -o output.avi mf://\*.png
