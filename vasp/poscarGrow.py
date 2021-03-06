#!/usr/bin/python

import sys,os
from scipy import *
from numpy import linalg
#mine
import poscarIO

def poscarGrow(poscarName,outputName,cx,cy,cz):
    poscar=open(poscarName,"r").readlines()
    nCopies=cx*cy*cz
 
    [basis,atypes,atoms,head,poscar] = poscarIO.read(poscar,frac_coord=True)
    natoms=len(atoms)
    v1,v2,v3=map(array,basis)
    ax,ay,az=map(array,zip(*atoms))

    #duplicate and tile
    ax=array(split(tile(ax,nCopies),nCopies))
    ay=array(split(tile(ay,nCopies),nCopies))
    az=array(split(tile(az,nCopies),nCopies))

    #Shift tiles in x dirn
    shift=array([1]*natoms)
    shiftfull=array([shift*j for j in range(cx)])
    for i in range(cy*cz):
        ax[i*cx+array(range(cx))]+=shiftfull

    #Shift tiles in y dirn
    shift=array([1]*natoms)
    shiftfull=array([shift*(j/cx) for j in range(cy*cx)])
    for i in range(cz):
        ay[i*cx*cy+array(range(cy*cx))]+=shiftfull

    #Shift tiles in z dirn
    shift=array([1]*natoms)
    shiftfull=array([shift*(j/cx/cy) for j in range(cy*cx*cz)])
    az+=shiftfull

    ax=ax.reshape([nCopies*natoms])
    ay=ay.reshape([nCopies*natoms])
    az=az.reshape([nCopies*natoms])
    
    #need to normalize with new positions with respect to duplicated lattice params
    v1*=cx
    v2*=cy
    v3*=cz
    basis=[v1,v2,v3]
    A=matrix([[float(cx),0.0,0.0],[0.0,float(cy),0.0],[0.0,0.0,float(cz)]])
    ax,ay,az=zip(*[linalg.solve(A,array(p).T)[:] for p in zip(ax,ay,az)])

#    if max(ax)>1.0 or max(ay)>1.0 or max(az)>1.0:
#        continue
#        exit(0)

    atoms=zip(ax,ay,az)

    #Sort atoms by type, as required by POSCAR
    types=list()
    j=0
    for i in atypes:
        types+=[j+1]*i
        j+=1
    newtypes=types*cx*cy*cz

    atoms,types=zip(*sorted(zip(atoms,newtypes),key=lambda x:x[1]))
    atypes=[i*nCopies for i in atypes]

    poscarIO.write(outputName,basis,atoms,atypes,head)
    return len(atoms),basis

if __name__=="__main__":

    def usage():
        print "Usage: %s <Input poscar/contcar> <Output> <dupl-x> <dupl-y> <dupl-z>"%sys.argv[0].split("/")[-1]
        print "Duplicates the poscar along the primary axes defined by the bounding box"

    if len(sys.argv)!=6:
        usage()
        exit(0)

    poscarName=sys.argv[1]
    outputName=sys.argv[2]

    cx,cy,cz=map(int,sys.argv[3:6])
    if cx<1 or cy<1 or cz<1:
        print "Cell duplication values must be >=1."
        exit(0)
    natoms,basis=poscarGrow(poscarName,outputName,cx,cy,cz)
    
    print "%d Atoms in grown cell"%(natoms)
    print "Bounds of simulation area are:"
    for i in basis:
        print i

