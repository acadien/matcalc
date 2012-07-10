#!/usr/bin/python

import sys
#mine
from poscarIO import readposcar
from orderParams import bondOrient
from numpy import *

poscar = open(sys.argv[1],"r").readlines()
[v1,v2,v3,atypes,ax,ay,az,head,poscar] = readposcar(poscar)
atoms=array(zip(ax,ay,az))
basis=array([v1,v2,v3])
bounds=array([v1[0],v2[1],v3[2]])

print bondOrient(atoms,basis,6)

