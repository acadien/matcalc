#!/usr/bin/python

import sys,os
from scipy import *
from numpy import linalg
#mine
import poscarIO

def poscarGrow(poscar,outputname,cx,cy,cz):
    nCopies=cx*cy*cz
 
    [v1,v2,v3,atypes,ax,ay,az,head,poscar] = poscarIO.read(poscar,frac_coord=True)
    natoms=len(ax)
    v1,v2,v3=map(array,[v1,v2,v3])
    ax,ay,az=map(array,[ax,ay,az])
    
    ax=array(split(tile(ax,nCopies),nCopies))
    ay=array(split(tile(ay,nCopies),nCopies))
    az=array(split(tile(az,nCopies),nCopies))

    #shift=array([v1[0]+v2[0]+v3[0]]*natoms)
    shift=array([1]*natoms)
    shiftfull=array([shift*j for j in range(cx)])
    for i in range(cy*cz):
        ax[i*cx+array(range(cx))]+=shiftfull

        #shift=array([v1[1]+v2[1]+v3[1]]*natoms)
    shift=array([1]*natoms)
    shiftfull=array([shift*(j/cx) for j in range(cy*cx)])
    for i in range(cz):
        ay[i*cx*cy+array(range(cy*cx))]+=shiftfull

        #shift=array([v1[2]+v2[2]+v3[2]]*natoms)
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
    atoms=zip(ax,ay,az)

    atypes=[i*nCopies for i in atypes]

    poscarIO.write(outputname,basis,atoms,atypes,head)
    return len(atoms),basis

if __name__=="__main__":

    def usage():
        print "Usage: %s <poscar/contcar> <dupl-x> <dupl-y> <dupl-z>"%sys.argv[0]
        print "Duplicates the poscar along the primary axes defined by the bounding box"

    if len(sys.argv)<4:
        usage()
        exit(0)

    outputname=sys.argv[1]+"_grown"
    poscar=open(sys.argv[1],"r").readlines()
    cx,cy,cz=map(int,sys.argv[2:5])
    if cx<1 or cy<1 or cz<1:
        print "Cell duplication values must be >=1."
        exit(0)
    natoms,basis=poscarGrow(poscar,outputname,cx,cy,cz)
    
    print "%d Atoms in grown cell"%(natoms)
    print "Bounds of simulation area are:"
    for i in basis:
        print i

