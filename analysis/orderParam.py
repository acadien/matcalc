#!/usr/bin/python

from scipy import special
import numpy as np
import math
from operator import mul
import pylab as pl
#mine
from struct_tools import sphang,neighborBasis,neighborOrtho
from rdf import rdf
from datatools import flatten

#Raw functions for calculating order parameters such as bond orientation, translational and tetrahedral ordering.

#Helper function for bond-orientation
#the spherical harmonics for the angle between two 3D points a and b
def pairSphereHarms(a,b,l):
    theta,phi=sphang(a,b)
    return special.sph_harm(np.array(range(-l,l+1)),l,theta,phi)
    
#bond-orientational: Q_l for each atom.  atomi>-1 selects a specific atom
def bondOrientation(atoms,bounds,l,ortho=True,atomi=-1,neighbCut=5.0,neighbs=None):
    #First make the neighbor list
    if neighbs==None:
        if ortho:
            neighbs=neighborOrtho(atoms,bounds,neighbCut,style="full")
        else:
            neighbs=neighborBasis(atoms,bounds,neighCut,style="full")
    print "starting Bond Orientation Calculation"
    #sum the spherical harmonic over ever neighbor pair
    if atomi==-1:
        Qlms = [sum( [ pairSphereHarms(atoms[i],atoms[j],l) for j in ineighbs ] ) \
                     for i,ineighbs in enumerate(neighbs) ] 
    #Grab the Qlm's for just 1 atom
    else:
        Qlms = [sum( [ pairSphereHarms(atoms[atomi],atoms[j],l) for j in neighbs[atomi] ] )]

    #Qlms[atom][m]

    print "before calc"
    Ql = [ (Qlm.conjugate()*Qlm *2*np.pi / (l+0.5)).sum()**0.5 for Qlm in Qlms] 
    print "done calc"
    return Ql

#translational: t
def translational(atoms,basis,neighbs=None,rcut=None):
    if rcut==None:
        rcut=basis[0][0]/2.0
    r,g = rdf(atoms,cutoff=rcut)
    h=map(math.fabs,g-1)
    
    tao = sum(h)/len(r)
    return tao
        
#tetrahedral
