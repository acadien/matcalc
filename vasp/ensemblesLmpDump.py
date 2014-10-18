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

possibleSuffix = [".rmsd",".cn",".tetra",".rmsd100",".rmsd1000",".rmsd5000",".rmsd10000",".2shelltetra",".2shellcn"]

utils.usage(["outcar"],1,1,"Note: Automatically looks for <OUTCAR.suffix> for per atom ensembles. to include in the dump")

outcarFilename = sys.argv[1]
dumpFilename = outcarFilename+"lmp.dump"
possibleSuffix = [i for i in possibleSuffix if os.path.isfile(outcarFilename+i)]
print possibleSuffix
ensembleHead   = " ".join([i.strip(".") for i in possibleSuffix])
ensembles = None
if len(possibleSuffix)>0:
    ensembles  = [parserGens.parseEnsemble(outcarFilename+i) for i in possibleSuffix]

nAtoms = outcarIO.nIons(outcarFilename)
basis = array(map(array,outcarIO.basis(outcarFilename)))
bounds = [[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
types = outcarIO.types(outcarFilename)

#Find the starting locations of atomic data in outcarfile
grepResults = subprocess.check_output("grep -b POSITION %s"%outcarFilename,shell=True).split("\n")
bytenums=[int(i.split(":")[0]) for i in grepResults if len(i)>2]
outcarParsed = parserGens.parseOutcarAtoms(bytenums,open(outcarFilename,"r"),nAtoms)

lammpsIO.writeDump(dumpFilename,basis,types,outcarParsed,ensembles,ensembleHead)


