#!/usr/bin/python

from scipy import special
from numpy import zeros,pi
#mine
from struct_tools import sphang,neighbors
from datatools import flatten

#Raw functions for calculating order parameters such as bond orientation, translational and tetrahedral ordering.

def pair_sphere_harm(atoma,atomb,l):
    theta,phi=sphang(atoma,atomb)
    return [special.sph_harm(m,l,theta,phi) for m in range(-l,l+1)]
            
#bond-orientational: Q_l
def bondOrient(atoms,basis,l,neighbs=None,rcut=5.0):

    #First make the neighbor list
    if neighbs==None:
        bounds=[[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
        neighbs=neighbors(atoms,bounds,rcut,style="half")

    #Average the spherical harmonic over ever neighbor pair
    Qlm = flatten([[pair_sphere_harm(atoms[i],atoms[j],l) for j in ineighbs] for i,ineighbs in enumerate(neighbs)])
    Qlmavg = map(lambda x: sum(x)/len(x),zip(*Qlm))
    Ql = 4.*pi*sum([abs(Qa.conjugate()*Qa) for Qa in Qlmavg])/(2*l+1.)
    return Ql

#translational: t
def translational(atoms,basis):
    pass
        
#tetrahedral
