#!/usr/bin/python

import sys,os
from scipy import *
#mine
import poscarIO
from duplicate import duplicate26

def poscar2qvoronoi(atoms,basis,atypes,qvfile="QV_input"):

    qvfp=open(qvfile,"w")

    NRealAtoms=len(atoms)

    #Ensure simulation bx is sufficiently orthorhombic
    if sum(map(fabs,basis[0]))-fabs(basis[0][0]) > 1e-10 or \
            sum(map(fabs,basis[1]))-fabs(basis[1][1]) > 1e-10 or \
            sum(map(fabs,basis[2]))-fabs(basis[2][2]) > 1e-10:
        print "Error: simulation axes are not sufficiently orthogonal.  Use script poscarRectify.py and regenerate doscar"
        exit(0)

    types=list()
    j=1
    for i in atypes:
        types+=[j]*i
        j+=1

    datoms,dtypes,dbasis=duplicate26(atoms,types,basis)

    bounds=[[-i/2.0,i*3.0/2] for i in [basis[0][0],basis[1][1],basis[2][2]]]
    def inBounds(point,bounds):
        for a,[b,c] in zip(point,bounds):
            if a>=c or a<b:
                return False
        return True

    datoms=[atom for atom in datoms if inBounds(atom,bounds)]
    Natoms=len(datoms)

    qvfp.write("3\n")
    qvfp.write("%d\n"%Natoms)
    for atom in datoms:
        qvfp.write("% 6.6f % 6.6f % 6.6f\n"%(atom[0],atom[1],atom[2]))
    qvfp.flush()
    qvfp.close()
    return qvfile,basis,datoms,NRealAtoms


if __name__=="__main__":

    def usage():
        print "Usage: %s <poscar/contcar> <QV filename>"%sys.argv[0]
        print "Makes a QV_input file from the poscar/contcar file.  Employs periodic boundary conditions."

    if len(sys.argv)<2:
        usage()
        exit(0)

    poscarFile=open(sys.argv[1],"r").readlines()
    basis,atypes,atoms,head,poscar = poscarIO.read(poscarFile)

    qvfile=None
    if len(sys.argv)>=3:
        qvfile=sys.argv[2]
        dummy,basis,datoms,nreal=poscar2qvoronoi(atoms,basis,atypes,qvfile)    
    else:
        dummy,basis,datoms,nreal=poscar2qvoronoi(atoms,basis,atypes)

    print "Bounds of simulation area are:"
    for i in basis:
        print i
    print "Make sure to remove excess atoms outside these bounds after analysis."
