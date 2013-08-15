#!/usr/bin/python

import sys
#mine
import poscarIO
import lammpsIO

def usage():
    print "%s <lammps dump> <config #> <POSCAR output>"%sys.argv[0].split("/")[-1]
    print "If only the dump file is provided, the number of configurations will be reported"
    print "Converts a lammps dump to POSCAR"

#If only the filename is provided, count the number of configurations and spit that out
if len(sys.argv)==2:
    dumpdata=open(sys.argv[1],"r").readlines()
    usage()
    print "\nNumber of configurations:",len([i for i in dumpdata if "ITEM: TIMESTEP" in i]),"\n"
    exit(0)    

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
poscarIO.write(posname,basis,atoms,atypes,header)
