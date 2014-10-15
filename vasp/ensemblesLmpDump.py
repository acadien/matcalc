#!/usr/bin/python

import sys
from math import *
import subprocess
import os
from numpy import array
#mine
import outcarIO
import lammpsIO

def parseEnsemble(ensembleFile):
    f = open(ensembleFile,"r")
    f.readline()
    for line in f:
        try: 
            yield map(float,line.split()[1:])
        except ValueError:
            pass

def parseOutcar(bytenums,outcarFile,nAtoms):
    for i,b in enumerate(bytenums):
        outcarFile.seek(b)
        outcarFile.readline()
        outcarFile.readline()
        atoms = [map(float,outcarFile.readline().split()[:3]) for a in range(nAtoms)]
        yield atoms

def usage():
    print "Usage: %s <outcar>"%sys.argv[0].split("/")[-1]
    print "Note: Automatically looks for <OUTCAR.suffix> for per atom ensembles. to include in the dump"

if len(sys.argv) not in [2]:
    usage()
    exit(0)

outcarFilename = sys.argv[1]
dumpFilename = outcarFilename+"lmp.dump"
possibleSuffix = [".rmsd",".cn",".tetra",".rmsd100",".rmsd1000",".rmsd5000",".rmsd10000"]
possibleSuffix = [i for i in possibleSuffix if os.path.isfile(outcarFilename+i)]
ensembleHead   = " ".join([i.strip(".") for i in possibleSuffix])
ensembles = None
if len(possibleSuffix)>0:
    ensembles  = [parseEnsemble(outcarFilename+i) for i in possibleSuffix]

nAtoms = outcarIO.nIons(outcarFilename)
basis = array(map(array,outcarIO.basis(outcarFilename)))
bounds = [[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
types = outcarIO.types(outcarFilename)

#Find the starting locations of atomic data in outcarfile
grepResults = subprocess.check_output("grep -b POSITION %s"%outcarFilename,shell=True).split("\n")
bytenums=[int(i.split(":")[0]) for i in grepResults if len(i)>2]
outcarParsed = parseOutcar(bytenums,open(outcarFilename,"r"),nAtoms)

lammpsIO.writeDump(dumpFilename,basis,types,outcarParsed,ensembles,ensembleHead)


