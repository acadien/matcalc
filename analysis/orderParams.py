#!/usr/bin/python

from scipy import special
from numpy import zeros,pi,nan
import math
from operator import mul
import pylab as pl
#mine
from struct_tools import sphang,neighbors
from paircor import paircor
from datatools import flatten

#Raw functions for calculating order parameters such as bond orientation, translational and tetrahedral ordering.

#bond-orientational: Q_l
def bondOrient(atoms,basis,l,neighbs=None,rcut=5.0):

    #Helper function for bond-orientation
    def pair_sphere_harm(atoma,atomb,l):
        theta,phi=sphang(atoma,atomb)
        return [special.sph_harm(m,l,theta,phi) for m in range(-l,l+1)]

    #First make the neighbor list
    if neighbs==None:
        bounds=[[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
        neighbs=neighbors(atoms,bounds,rcut,style="half")

    #Average the spherical harmonic over ever neighbor pair
    Qlm = flatten([[pair_sphere_harm(atoms[i],atoms[j],l) for j in ineighbs] for i,ineighbs in enumerate(neighbs)])
    Qlmavg = map(lambda x: sum(x)/len(x),zip(*Qlm))

    #Sum over m vals and sqrt.
    Ql = (4.*pi*sum([abs(Qa.conjugate()*Qa) for Qa in Qlmavg])/(2*l+1.))**0.5

    return Ql

#translational: t
def translational(atoms,basis,neighbs=None,rcut=None):
    if rcut==None:
        rcut=basis[0][0]/2.0
    r,g = paircor(atoms,cutoff=rcut)
    h=map(math.fabs,g-1)
    
    tao = sum(h)/len(r)
    return tao
        
#tetrahedral
