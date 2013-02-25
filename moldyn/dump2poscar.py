#!/usr/bin/python

import sys
#mine
import poscarIO
import lammpsIO

def usage():
    print "%s <lammps dump> <config #> <POSCAR output>"%sys.argv[0].split("/")[-1]
    print "Converts a lammps dump to POSCAR"

if len(sys.argv)!=4:
    usage()
    exit(0)

dumpdata=open(sys.argv[1],"r").readlines()
configN=int(sys.argv[2])
posname=sys.argv[3]

#Read in lammps data
basis,types,atoms,header=lammpsIO.dumpReadConfig(dumpdata,configN)

#Write out POSCAR data
atypes=[types.count(j) for j in sorted(list(set(types)))]
poscarIO.write(posname,basis,atoms,atypes,header,frac=True)
