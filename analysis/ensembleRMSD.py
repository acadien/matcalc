#!/usr/bin/python

import sys
from math import *
import subprocess
import os
from numpy import *
#mine
import outcarIO
import lammpsIO
import rootMeanSquareDist
import utils

utils.usage(["<Outcar/Lammpsdump>"],1,1)

filename = sys.argv[1]
header=["AverageRMSD PerAtomRMSD\n"]
rmsddata = header

outcarFlag=False
lammpsFlag=False
if "OUTCAR" in filename:
    outcarFlag=True
    lammpsFlag=False
else:
    outcarFlag=False
    lammpsFlag=True

atoms=list()
if outcarFlag:
    nAtoms = outcarIO.nAtoms(filename)
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
delT,rmsd,rmsdByAtom=rootMeanSquareDist.rootMeanSquareDistRef(atoms,0,basis,byAtom=True)
rmsddata+=[str(y)+" "+" ".join(map(str,z))+"\n" for y,z in zip(rmsd,rmsdByAtom)]

rmsdfile=filename+".rmsd"

print "Writing %s."%rmsdfile
open(rmsdfile,"w").writelines(rmsddata)
