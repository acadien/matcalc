#!/usr/bin/python

from scipy import special
import numpy as np
import math
from operator import mul
import pylab as pl
#mine
from voronoiNeighbors import voronoiNeighbors
from struct_tools import sphang,neighborBasis,neighborOrtho,dist_periodic
from rdf import rdf
from datatools import flatten,windowAvg

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


def coordinationNumbers(atoms,basis,ortho=True,neighbs=None,rcut=None):
    if rcut==None:
        #Set rcut to be the first minimum of g(r)
        rvals,gr = rdf(atoms,cutoff=6.0)

        #Smoothed G(r)
        sgr=windowAvg(gr,n=100)
        #derivative of smoothed-G(r)
        dsgr = windowAvg([sgr[i+1]-sgr[i] for i in range(len(sgr)-1)],n=50) 
        
        #Find the first minima by chopping at the first negative and the 2nd positive.
        first_neg = [i for i,v in enumerate(dsgr) if v<0][0]
        dsgr = dsgr[first_neg:]
        second_pos = [i for i,v in enumerate(dsgr) if v>0][0]
        dsgr = dsgr[:second_pos]
        rcut = rvals[dsgr.index(min(dsgr))+first_neg ]
        print rcut
#        rcut=rcut*2
#        pl.plot(rvals,gr)
#        pl.plot(rvals,[i for i in sgr],lw=3)
#        pl.plot([rcut,rcut],[min(sgr),max(sgr)])
#        pl.show()

    #First make the neighbor list
    if neighbs==None:
#        kwargs={"atoms":atoms,"basis":basis,"atypes":[len(atoms)],"style":"full"}
#        neighbs=voronoiNeighbors(atoms=atoms,basis=basis,atypes=[len(atoms)],style="full")
        bounds=[[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
#        neighbs=neighborOrtho(atoms,bounds,rcut,style="full")
        neighbs=neighborBasis(atoms,basis,rcut,style="full")
    cns = [sum([1 for j in ineighbs if dist_periodic(atoms[i],atoms[j],np.array(zip(*bounds)[1])) < rcut])\
                    for i,ineighbs in enumerate(neighbs)]

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
