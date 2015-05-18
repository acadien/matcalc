#!/usr/bin/python

import sys
from math import *
import subprocess
import os
from scipy import array,zeros
import pylab as pl
#mine
import outcarIO
import lammpsIO
import orderParam
import utils 

utils.usage(["<outcar> <l> -2neighb"],2,3,"-2neighb: 2nd shell neighbors are used when calculating Ql")

rcut = None
filename = sys.argv[1]
lval = int(sys.argv[2])
sh2 = False
if "-2neighb" in sys.argv:
    sh2=True

outcarFlag=False
lammpsFlag=False
if "OUTCAR" in filename:
    outcarFlag=True
    lammpsFlag=False
else:
    outcarFlag=False
    lammpsFlag=True


if outcarFlag:
    nAtoms = outcarIO.nAtoms(filename)
    basis = array(map(array,outcarIO.basis(filename)))

    #Find the starting locations of atomic data in outcarfile
    grepResults = subprocess.check_output("grep -b POSITION %s"%filename,shell=True).split("\n")
    bytenums=[int(i.split(":")[0]) for i in grepResults if len(i)>2]

    outcar= open(filename,"r")
    bondOrderOut =open(filename+".BO%d"%lval,"w")
    if sh2:
        bondOrderOut =open(filename+".BO%dsh2"%lval,"w")

    bondOrderOut.write("AverageBO%d PerAtomBO%d\n"%(lval,lval))
    for i,b in enumerate(bytenums):
        outcar.seek(b)
        outcar.readline()
        outcar.readline()

        atoms = [map(float,outcar.readline().split()[:3]) for a in range(nAtoms)]
        if sh2:
            q,rcut = orderParam.bondOrientation2sh(atoms,basis,lval)
        else:
            q,rcut = orderParam.bondOrientation(atoms,basis,lval)

        qAvg = sum(q)/len(q)
        qs = " ".join(map(str,q))
        bondOrderOut.write(str(qAvg)+" "+qs+"\n")
        #print i

if lammpsFlag:
    nAtoms = lammpsIO.nAtoms(filename)
    basisByteNums = lammpsIO.basisBytes(filename)
    atomsByteNums = lammpsIO.atomsBytes(filename)
    
    bondOrderOut = open(filename+".BO%d"%lval,"w")
    if sh2:
        bondOrderOut =open(filename+".BO%dsh2"%lval,"w")

    bondOrderOut.write("AverageBO%d PerAtomBO%d\n"%(lval,lval))
    for i,(bByte,aByte) in enumerate(zip(basisByteNums,atomsByteNums)):
        basis = lammpsIO.parseBasis(filename,bByte)

        atoms,dummy = lammpsIO.parseAtoms(filename,aByte,nAtoms,basis)

        if sh2:
            q,rcut = orderParam.bondOrientation2sh(atoms,basis,lval)
        else:
            q,rcut = orderParam.bondOrientation(atoms,basis,lval)

        qAvg = sum(q)/len(q)
        qs = " ".join(map(str,q))
        bondOrderOut.write(str(qAvg)+" "+qs+"\n")
        #print i
