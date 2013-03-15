#!/usr/bin/python

from scipy import special
import numpy as np
import math
from operator import mul
import pylab as pl
#mine
from voronoiNeighbors import voronoiNeighbors
from struct_tools import *
from rdf import rdf
from datatools import flatten,windowAvg

#Raw functions for calculating order parameters such as bond orientation, translational and tetrahedral ordering.

#Helper function for bond-orientation
#the spherical harmonics for the angle between two 3D points a and b
def pairSphereHarms(a,b,l):
    theta,phi=sphang(a,b)
    return special.sph_harm(np.array(range(-l,l+1)),l,theta,phi)

#Given a set of atoms, finds the first minimum of G(r) and 
#returns neighbors with bonds shorter than that length
def neighborHelperPeriodic(atoms,basis,rcut=None):
    natoms = len(atoms)
    allatoms = makeGhosts(atoms,basis)

    if rcut==None:
        #Set rcut to be the first minimum of g(r)
        rvals,gr = rdf(atoms,cutoff=6.0)

        #Smoothed G(r)
        sgr=windowAvg(gr,n=25)
        #derivative of smoothed-G(r)
        dsgr = windowAvg([sgr[i+1]-sgr[i] for i in range(len(sgr)-1)],n=50) 
        
        #Find the first minima by chopping at the first negative and the 2nd positive.
        first_neg = [i for i,v in enumerate(dsgr) if v<0][0]
        dsgr = dsgr[first_neg:]
        second_pos = [i for i,v in enumerate(dsgr) if v>0][0]
        dsgr = dsgr[:second_pos]
        rcut = rvals[dsgr.index(min(dsgr))+first_neg ]

    #Should probably check out this plot before continuing
    #print rcut
    #pl.plot(rvals,gr)
    #pl.plot(rvals,[i for i in sgr],lw=3)
    #pl.plot([rcut,rcut],[min(sgr),max(sgr)])
    #pl.show()
    neighbs = neighborBasis(allatoms,natoms,basis,rcut)
    return neighbs,allatoms


#bond-orientational: Q_l for each atom.  atomi>-1 selects a specific atom
def bondOrientation(atoms,basis,l,neighbs=None,rcut=None):
    #First make the neighbor list
    if neighbs==None:
        neighbs,allatoms = neighborHelperPeriodic(atoms,basis,rcut)

    #sum the spherical harmonic over ever neighbor pair
    Qlms = [sum( [ pairSphereHarms(atoms[i],allatoms[j],l) for j in ineighbs ] ) \
                for i,ineighbs in enumerate(neighbs) ] 

    Ql = [ ((Qlm.conjugate()*Qlm *2*np.pi / (l+0.5)).sum()**0.5).real for Qlm in Qlms] 

    return Ql


def coordinationNumber(atoms,basis,l=None,neighbs=None,rcut=None):
    #l: not used
    
    if neighbs==None:
        neighbs,allatoms = neighborHelperPeriodic(atoms,basis,rcut)

    cns = map(len,neighbs)
        
    return cns
            

#translational: t
def translational(atoms,basis,neighbs=None,rcut=None):
    if rcut==None:
        rcut=basis[0][0]/2.0
    r,g = rdf(atoms,cutoff=rcut)
    h=map(math.fabs,g-1)
    
    tao = sum(h)/len(r)
    return tao
        
#tetrahedral
