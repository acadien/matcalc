#!/usr/bin/python
import sys
from math import *
#mine
import poscarIO

#Reading in the options and preparing for file reading
if len(sys.argv)<3:
    print "Usage:"
    print sys.argv[0]+" <file with a bunch of poscars> <suffix> <Atom type1> <Atom type2> ..."
    exit(0)

poscar = open(sys.argv[1],"r").readlines()

suffix = sys.argv[2]

atoms=" ".join(sys.argv[3:])+"   "

ind=0
while True:
    if len(poscar)<=1:
        break

    poscar,poscarout = poscarIO.split(poscar)
    poscarout[0]=atoms+poscarout[0]
    if poscarout==poscar==-1:
        break

    fname="POSCAR"+suffix+"%3.3d"%ind
    open(fname,"w").writelines(poscarout)
    ind+=1
