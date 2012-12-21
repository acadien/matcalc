#!/usr/bin/python

import sys
from numpy import *
#mine
import poscarIO
from orderParams import bondOrient,translational

def usage():
    print "%s <POSCAR> <l for bond orientation>"%sys.argv[0].split("/")[-1]

if len(sys.argv)<2:
    usage()
    exit(0)

poscar = open(sys.argv[1],"r").readlines()
[basis,atypes,atoms,head,poscar] = poscarIO.read(poscar)
bounds=array([basis[0][0],basis[1][1],basis[2][2])

print "Orientational Param:", bondOrient(atoms,basis,int(sys.argv[2]))
print "Translational Param:", translational(atoms,basis)
