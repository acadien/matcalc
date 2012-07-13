#!/usr/bin/python

import sys
#mine
from poscarIO import readposcar
from orderParams import bondOrient,translational
from numpy import *

def usage():
    print "%s <POSCAR> <l for bond orientation>"%sys.argv[0].split("/")[-1]

if len(sys.argv)<2:
    usage()
    exit(0)

poscar = open(sys.argv[1],"r").readlines()
[v1,v2,v3,atypes,ax,ay,az,head,poscar] = readposcar(poscar)
atoms=array(zip(ax,ay,az))
basis=array([v1,v2,v3])
bounds=array([v1[0],v2[1],v3[2]])

print "Orientational Param:", bondOrient(atoms,basis,int(sys.argv[2]))
print "Translational Param:", translational(atoms,basis)
