#!/usr/bin/python

from scipy import special
import numpy as np
import math
from operator import mul
import pylab as pl
#mine
from neighbors import neighbors
from struct_tools import *
from rdf import rdf
from datatools import flatten,windowAvg

#Raw functions for calculating order parameters such as bond orientation, translational and tetrahedral ordering.

#Helper function for bond-orientation
#the spherical harmonics for the angle between two 3D points a and b
def pairSphereHarms(a,b,l):
    theta,phi=sphang(a,b)
    return special.sph_harm(np.array(range(-l,l+1)),l,theta,phi)

#bond-orientational: Q_l for each atom.  atomi>-1 selects a specific atom
def bondOrientation(atoms,basis,l,neighbs=None,rcut=None):
    #First make the neighbor list
    if neighbs==None:
        if rcut==None:
            rcut = generateRCut(atoms,debug=True)
        bounds=[[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
        neighbs = neighbors(atoms,bounds,rcut)
#        allatoms,reduceGhost = makeGhosts(atoms,basis)

    #sum the spherical harmonic over ever neighbor pair
    Qlms = [sum( [ pairSphereHarms(atoms[i],minImageAtom(atoms[i],atoms[j],basis),l) for j in ineighbs ] ) \
                for i,ineighbs in enumerate(neighbs) ] 
    Ql = [ ((Qlm.conjugate()*Qlm *2*np.pi / (l+0.5)).sum()**0.5).real for Qlm in Qlms] 

    return Ql


def coordinationNumber(atoms,basis,l=None,neighbs=None,rcut=None):
    #l: not used
    
    if neighbs==None:
        if rcut==None:
            rcut = generateRCut(atoms,debug=True)
        bounds=[[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
        neighbs = neighbors(atoms,bounds,rcut)

    cns = map(len,neighbs)
        
    return cns
            

#translational: t
def translational(atoms,basis,neighbs=None,rcut=None):
    if rcut==None:
        rcut = generateRCut(atoms,debug=True)
    r,g = rdf(atoms,cutoff=rcut)
    h=map(math.fabs,g-1)
    
    tao = sum(h)/len(r)
    return tao
        
#tetrahedral
