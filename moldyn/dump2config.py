#!/usr/bin/python

import sys
#mine
from lammpsIO import *

def usage():
    print sys.argv[0].split("/")[-1]+" <input dump file> <configuration #> <output config file>"
    print "LAMMPS dump file must have bounds, atom type and atomic locations."
    

if len(sys.argv)!=4:
    usage()
    exit(0)

iDump=open(sys.argv[1],"r").readlines()
cfg=int(sys.argv[2])
bounds,types,atoms,head=dumpReadConfig(iDump,cfg)
oCnfg=open(sys.argv[3],"w")
dumpWriteConfig(oCnfg,bounds,types,atoms,head)

