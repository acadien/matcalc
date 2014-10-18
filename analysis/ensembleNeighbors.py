#!/usr/bin/python

import sys
from math import *
import subprocess
import os
from numpy import array
#mine
import outcarIO,xdatcarIO,poscarIO
import lammpsIO
import neighbors

def usage():
    print "Usage: %s <Outcar/Lammpsdump/Xdatcar> <rcut> <poscarfile if xdat used>"%sys.argv[0].split("/")[-1]

if len(sys.argv) not in [3,4]:
    usage()
    exit(0)

filename = sys.argv[1]
rcut = float(sys.argv[2])
outcarFlag=False
lammpsFlag=False
xdatcarFlag=False
if "OUTCAR" in filename:
    outcarFlag=True
elif "XDATCAR" in filename:
    xdatcarFlag=True
    poscarFile = sys.argv[3]
else:
    lammpsFlag=True

neighbsfile=filename+".neighb"
header=["Spaces Seperate Neighbs, Commas Seperate Atoms, Lines Seperate Arrangements\n"]
lines=header

if xdatcarFlag:
    basis,dummy,atoms,dummy,dummy = poscarIO.read(poscarFile)
    nAtoms = len(atoms)
    bounds = [[0,1],[0,1],[0,1]]
    rcut /= basis[0][0]
    for i,atoms in enumerate(xdatcarIO.read2(filename)):
        neighbs = neighbors.neighbors(atoms,bounds,rcut)
        lines += [",".join([" ".join(map(str,atomn)) for atomn in neighbs])+"\n"]

if outcarFlag:
    nAtoms = outcarIO.nIons(filename)
    basis = array(map(array,outcarIO.basis(filename)))
    bounds = [[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]

    #Find the starting locations of atomic data in outcarfile
    grepResults = subprocess.check_output("grep -b POSITION %s"%filename,shell=True).split("\n")
    bytenums=[int(i.split(":")[0]) for i in grepResults if len(i)>2]

    outcar = open(filename,"r")
    for i,b in enumerate(bytenums):
        outcar.seek(b)
        outcar.readline()
        outcar.readline()
        atoms = [map(float,outcar.readline().split()[:3]) for a in range(nAtoms)]
        neighbs = neighbors.neighbors(atoms,bounds,rcut)
        lines += [",".join([" ".join(map(str,atomn)) for atomn in neighbs])+"\n"]

if lammpsFlag:
    nAtoms = lammpsIO.nAtoms(filename)
    basisByteNums = lammpsIO.basisBytes(filename)
    atomsByteNums = lammpsIO.atomsBytes(filename)
    
    for i,(bByte,aByte) in enumerate(zip(basisByteNums,atomsByteNums)):
        basis = lammpsIO.parseBasis(filename,bByte)
        bounds = [[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
        atoms,dummy = lammpsIO.parseAtoms(filename,aByte,nAtoms,basis)
        neighbs = neighbors.neighbors(atoms,bounds,rcut)
        lines += [",".join([" ".join(map(str,atomn)) for atomn in neighbs])+"\n"]

print "Writing %s."%neighbsfile
open(neighbsfile,"w").writelines(lines)
