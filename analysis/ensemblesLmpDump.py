#!/usr/bin/python

import sys
from math import *
import subprocess
import os
from numpy import array
#mine
import outcarIO
import lammpsIO
import parserGens
import utils

possibleSuffix = [".rmsd",".cn",".tetra",".rmsd100",".rmsd1000",".rmsd5000",".rmsd10000",".2shelltetra",".2shellcn",".2shellrmsd100"]

utils.usage(["<outcar,lamps dump>"],1,1,"Note: Automatically looks for <OUTCAR/dump.suffix> for per atom ensembles. to include in the dump")

fname = sys.argv[1]
dumpFilename = fname+"lmp.dump"
possibleSuffix = [i for i in possibleSuffix if os.path.isfile(fname+i)]
print possibleSuffix
ensembleHead   = " ".join([i.strip(".") for i in possibleSuffix])
ensembles = None
if len(possibleSuffix)>0:
    ensembles  = [parserGens.parseEnsemble(fname+i) for i in possibleSuffix]

if "OUTCAR" in fname:
    nAtoms = outcarIO.nAtoms(fname)
    basis = array(map(array,outcarIO.basis(fname)))
    bounds = [[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
    types = outcarIO.types(fname)

    #Find the starting locations of atomic data in outcarfile
    grepResults = subprocess.check_output("grep -b POSITION %s"%fname,shell=True).split("\n")
    bytenums=[int(i.split(":")[0]) for i in grepResults if len(i)>2]
    outcarParsed = parserGens.parseOutcarAtoms(bytenums,open(fname,"r"),nAtoms)

    lammpsIO.writeDump(dumpFilename,basis,types,outcarParsed,ensembles,ensembleHead)

else:
    nAtoms = lammpsIO.nAtoms(fname)
    basis = lammpsIO.basis(fname)
    bounds = [[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
    b = lammpsIO.atomsBytes(fname)
    dummy,types = lammpsIO.parseAtoms(fname,b[0],nAtoms,basis)
    dumpParsed = parserGens.parseLammpsAtoms(b,fname,nAtoms)

    lammpsIO.writeDump(dumpFilename,basis,types,dumpParsed,ensembles,ensembleHead)
