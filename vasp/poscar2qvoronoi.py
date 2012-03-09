#!/usr/bin/python

import sys,os
from scipy import *
#mine
from poscarIO import readposcar
from duplicate import duplicate26

def poscar2qvoronoi(**kwargs):
    #Possible args: 
    #Option1: poscar, [qvfile]
    #Option2: atoms, basis, atypes, [qvfile]
    if "poscar" in kwargs:
        poscar=kwargs["poscar"]
        [v1,v2,v3,atypes,ax,ay,az,head,poscar] = readposcar(poscar)        
        basis=[v1,v2,v3]
        atoms=zip(ax,ay,az)
    else:
        atypes=kwargs["atypes"]
        basis=kwargs["basis"]
        [v1,v2,v3]=basis
        atoms=kwargs["atoms"]
    if "qvfile" in kwargs:
        qvfile=kwargs["qvfile"]
    else:
        qvfile="QV_input"

    qvfp=open(qvfile,"w")

    NRealAtoms=len(atoms)

    #Ensure simulation bx is sufficiently orthorhombic
    if sum(map(fabs,v1))-fabs(v1[0]) > 1e-10 or \
            sum(map(fabs,v2))-fabs(v2[1]) > 1e-10 or \
            sum(map(fabs,v3))-fabs(v3[2]) > 1e-10:
        print "Error: simulation axes are not sufficiently orthogonal.  Use script poscarRectify.py and regenerate doscar"
        exit(0)

    types=list()
    j=1
    for i in atypes:
        types+=[j]*i
        j+=1

    datoms,dtypes,dbasis=duplicate26(atoms,types,basis)

    bounds=[[-i/2.0,i*3.0/2] for i in [v1[0],v2[1],v3[2]]]
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

    poscar=open(sys.argv[1],"r").readlines()

    qvfile=None
    if len(sys.argv)>=3:
        qvfile=sys.argv[2]
    
    if qvfile==None:
        dummy,basis,datoms,nreal=poscar2qvoronoi(poscar=poscar)
    else:
        dummy,basis,datoms,nreal=poscar2qvoronoi(poscar,qvfile)

    print "Bounds of simulation area are:"
    for i in basis:
        print i
    print "Make sure to remove excess atoms outside these bounds after analysis."
