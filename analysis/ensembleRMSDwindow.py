#!/usr/bin/python

import sys
from math import *
import subprocess
import os
#mine
import utils
import outcarIO
import lammpsIO
import rootMeanSquareDist

from scipy import array,zeros
import pylab as pl

utils.usage(["<Outcar/Lammpsdump>","<windowsz>"],2,2)

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

window = int(sys.argv[2])

atoms=list()
if outcarFlag:
    nAtoms = outcarIO.nIons(filename)
    basis = array(map(array,outcarIO.basis(filename)))

    #Find the starting locations of atomic data in outcarfile
    grepResults = subprocess.check_output("grep -b POSITION %s"%filename,shell=True).split("\n")
    bytenums=[int(i.split(":")[0]) for i in grepResults if len(i)>2]

    outcar = open(filename,"r")
    for i,b in enumerate(bytenums):
        outcar.seek(b)
        outcar.readline()
        outcar.readline()
        atoms.append([map(float,outcar.readline().split()[:3]) for a in range(nAtoms)])

if lammpsFlag:
    nAtoms = lammpsIO.nAtoms(filename)
    basisByteNums = lammpsIO.basisBytes(filename)
    atomsByteNums = lammpsIO.atomsBytes(filename)
    
    for i,(bByte,aByte) in enumerate(zip(basisByteNums,atomsByteNums)):
        basis = lammpsIO.parseBasis(filename,bByte)
        try:
            a,dummy = lammpsIO.parseAtoms(filename,aByte,nAtoms,basis)
        except IndexError:
            break
        atoms.append(a)

atoms = array(atoms)
delT,rmsd,rmsdByAtom=rootMeanSquareDist.rootMeanSquareDistWindow(atoms,window,basis,byAtom=True)

rmsdfile=filename+".rmsd%d"%window
header=["AverageRMSD PerAtomRMSD\n"]
rmsddata=header+[str(y)+" "+" ".join(map(str,z))+"\n" for y,z in zip(rmsd,rmsdByAtom)]

print "Writing %s."%rmsdfile
open(rmsdfile,"w").writelines(rmsddata)
