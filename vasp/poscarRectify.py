#!/usr/bin/python

import pylab as pl
import sys
from scipy import *
#mine
from poscarsPlot import plotsimulation
from poscarIO import *
from struct_tools import *

def rectify_basis(basis,atoms):
    #Given: 3 vectors defining the boundaries of a simulation
    #And a collection of atoms that are fractional combinations of basis
    #Returns: Basis orthognoal vectors such with same atoms 

    xprimary=array([1.,0.,0.])
    yprimary=array([0.,1.,0.])
    zprimary=array([0.,0.,1.])
    print basis
    [a,b,c]=basis

    #Rotate the whole system so that a lies in [1,0,0] direction
    M1=rotmatx(a,xprimary)
    a=dot(M1,a)
    b=dot(M1,b)
    c=dot(M1,c)
    atoms=[dot(M1,atom) for atom in atoms]

    """
    #Trim off x part of b & c
    for atom in atoms:
        if atom[0]>a[0]:
            atom[0]-=a[0]
        if atom[0]<0:
            atom[0]+=a[0]
    b[0]=0
    c[0]=0

    #Rotate the whole system so that a remains in [1,0,0] and b lies in [0,1,0]
    M2=rotmatx(b,yprimary)
    a=dot(M2,a)
    b=dot(M2,b)
    c=dot(M2,c)
    atoms=[dot(M2,atom) for atom in atoms]

    #Trim off y part of c
    for atom in atoms:
        if atom[1]>b[1]:
            atom[1]-=b[1]
        if atom[1]<0:
            atom[1]+=b[1]
    c[1]=0
    """
    basis=[a,b,c]    

    return basis,atoms

def plotbasis(basis,aa=None,c='r',ls='-',shift=None):

    [v1,v2,v3]=basis

    if shift==None:
        shift=zeros([3,3])

    if aa==None:
        fig=pl.figure()
        aa = fig.add_subplot(111,projection='3d')

    end=v1+v2+v3
    v12=v1+v2
    v13=v1+v3
    v23=v2+v3
    for vec in [v1,v2,v3]:
        aa.plot(*zip(shift+[0,0,0],shift+vec),c=c,ls=ls)
    for vec in [v13,v23,v12]:
        aa.plot(*zip(shift+vec,shift+end),c=c,ls=ls)
    aa.plot(*zip(shift+v1,shift+v12),c=c,ls=ls)
    aa.plot(*zip(shift+v1,shift+v13),c=c,ls=ls)
    aa.plot(*zip(shift+v2,shift+v12),c=c,ls=ls)
    aa.plot(*zip(shift+v2,shift+v23),c=c,ls=ls)
    aa.plot(*zip(shift+v3,shift+v13),c=c,ls=ls)
    aa.plot(*zip(shift+v3,shift+v23),c=c,ls=ls)

def plotatoms(atoms,types,aa=None):

    if aa==None:
        fig=pl.figure()
        aa = fig.gca(projection='3d')

    colors=["red","blue","green","yellow","orange","purple","black"]

    ax,ay,az=zip(*atoms)
    ind=0
    for cindex,i in enumerate(types):
        aa.scatter(ax[ind:ind+i],ay[ind:ind+i],az[ind:ind+i],c=colors[cindex],marker="o")
        ind+=i
    return aa

def plotvectors(L):
    fig=pl.figure()
    aa = fig.add_subplot(111,projection='3d')
    for vec in L:
        aa.plot([0,vec[0]],[0,vec[1]],[0,vec[2]])

if __name__=="__main__":

    if len(sys.argv)<2:
        print "Usage:"
        print sys.argv[0]+" <poscar> <1:overwrite poscar>"
        exit(0)    

    poscar = open(sys.argv[1],"r").readlines() 
    writeenabled=False
    if len(sys.argv) == 3:
        if sys.argv[2]=="1":
            writeenabled = True 

    [v1,v2,v3,atypes,ax,ay,az,head,poscar] = readposcar(poscar)
    
    basis=[v1,v2,v3]
    rbasis=map(array,basis)
    
    atoms=zip(ax,ay,az)
    ratoms=map(array,atoms)

    rbasis,ratoms=rectify_basis(rbasis,ratoms)

    if writeenabled:
        val=raw_input("\nHit any-key to over-write POSCAR, ctrl+d to cancel.")
        writeposcar(sys.argv[1],rbasis,ratoms,atypes,head)

    fig1=pl.figure()
    aa1=plotsimulation(basis,atoms,atypes,fig1)
    pl.title("Original")
    fig2=pl.figure()
    aa2=plotsimulation(rbasis,ratoms,atypes,fig2)
    pl.title("Modified")
    pl.show()

