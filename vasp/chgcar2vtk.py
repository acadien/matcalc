#!/usr/bin/python

import sys
#mine
from chgcarIO import *

def usage():
    print "Usage: %s <input chgcar> <optional:output vtk filename>"%sys.argv[0].split("/")[-1]

if len(sys.argv)<2:
    usage()
    exit(0)
    
chgfname=sys.argv[1]
if len(sys.argv)==3:
    vtkfname=sys.argv[2]
else:
    vtkfname=chgfname+".vtk"


#chgcar=open(chgfname).readlines()
poscardata,gridsz,chgdata=readchgcar(open(chgfname,"r").readlines())
writeVTK(vtkfname,poscardata,gridsz,chgdata,clean=True)
