#!/usr/bin/python

import sys
from numpy import *
#mine
import poscarIO
from lammpsIO import bounds2lohi
from struct_tools import flatten

def poscar2dump(pcarFile,lmpFile,scale=None):
    try:
        poscar=open(pcarFile,"r").readlines()
        lammps=open(lmpFile,"w")
    except IOError:
        print "Error: %s Unable to open POSCAR/Lammps file."%sys.argv[0].split("/")[-1]
        usage()
        exit(0)

    if scale==None:
        scale=1.0

    [basis,atypes,atoms,head,poscar] = poscarIO.read(poscar)

    basis*=scale

    #Convert from POSCAR style basis vectors to LAMMPS style boundaries.
    xhi,yhi,zhi,xy,xz,yz=bounds2lohi(basis)
    basis=matrix([[xhi,0,0],[xy,yhi,0],[xz,yz,zhi]])

    N=sum(atypes)
    Ntypes=len(atypes)
    types=[i for i in flatten([[i]*v for i,v in enumerate(atypes)])]
    data="#Generated by %s from %s\n\n"%(sys.argv[0].split("/")[-1],pcarFile)
    data+="%d atoms\n"%N
    data+="%d atom types\n\n"%Ntypes
    data+="0.0000  %4f  xlo xhi\n"%xhi
    data+="0.0000  %4f  ylo yhi\n"%yhi
    data+="0.0000  %4f  zlo zhi\n"%zhi
    data+="% 4f % 4f % 4f xy xz yz\n\n"%(xy,xz,yz)
    data+="Atoms\n\n"
    for i,(t,atom) in enumerate(zip(types,atoms)):
        x,y,z=atom#(atom*basis).tolist()[0]
        data+="\t %d \t %d \t % 5f\t% 5f\t% 5f\n"%(i+1,t+1,x,y,z)
    lammps.write(data)
    lammps.close()

def usage():
    print "%s <in:POSCAR-file> <out:LAMMPS config-file> <scaling factor A>"%(sys.argv[0].split("/")[-1])


if __name__=="__main__":
    if len(sys.argv)<3:
        usage()
        exit(0)
    
    if len(sys.argv)==3:
        poscar2dump(sys.argv[1],sys.argv[2])
    elif len(sys.argv)==4:
        poscar2dump(sys.argv[1],sys.argv[2],float(sys.argv[3]))
