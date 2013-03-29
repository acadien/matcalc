#!/usr/bin/python
#Raw functions for calculating order parameters such as bond orientation, translational and tetrahedral ordering.

import numpy as np
from scipy.special import sph_harm as sph_harm
from scipy import conj
import math
from operator import mul
import pylab as pl
#mine
from neighbors import neighbors,full2half
from struct_tools import *
from rdf import rdf,generateRCut,rdf_periodic

#local bond-orientational: Q_l for each atom.  atomi>-1 selects a specific atom
def bondOrientation(atoms,basis,l,neighbs=None,rcut=None,debug=True):

    if neighbs==None:
        if rcut==None:
            rcut = generateRCut(atoms,debug=debug)
        bounds=[[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
        neighbs = neighbors(atoms,bounds,rcut)

    #sum the spherical harmonic over ever neighbor pair
    Qlms = [sum( [ pairSphereHarms(atoms[i],minImageAtom(atoms[i],atoms[j],basis),l) for j in ineighbs ] ) \
                for i,ineighbs in enumerate(neighbs) ] 
    Ql = [ ((Qlm.conjugate()*Qlm *2*np.pi / (l+0.5)).sum()**0.5).real for Qlm in Qlms] 

    return Ql

#Helper function, returns Qlm values m=(-l .. 0 .. +l) for a specific atom pair: atomi,atomj
def bondOrientR(atoms,basis,l,atomi,atomj):
    ai = atoms[atomi]
    aj = atoms[atomj]
    Qlm = pairSphereHarms(ai,minImageAtom(ai,aj,basis),l)

    return Qlm

#Bond antle correlation function Gl as defined in:
#       Nature Materials, Vol2, Nov. 2003, Sastry & Angell
def bondAngleCorr(atoms,basis,l,neighbs=None,rcut=None,debug=False):
    
    print "Start Bond Angle Correlation Calculation"

    if neighbs==None:
        #if rcut==None:
        #    rcut = generateRCut(atoms,debug=debug)
        rcut = 6.0
        bounds=[[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
        #Half neighbor list
        hneighbs = neighbors(atoms,bounds,rcut,style="half")    

    #At distances rbins calculate the bond angle correlation function
    nbins = 256
    delr  = rcut/nbins

    #Histogram of bond lengths
    rbins = [i*delr for i in range(nbins)]
    bcnts = [0 for i in range(nbins)]
    gvals = [0.0 for i in range(nbins)]

    #Get the atomic pairs at each bond length
    for i,ineighbs in enumerate(hneighbs):
    #    print i,len(ineighbs)
        for j in ineighbs:
            #i & j make an atom pair, d is the bond length between them
            jatom = minImageAtom(atoms[i],atoms[j],basis)
            d = dist(atoms[i],jatom)
            bbin=int(d/delr)
            bcnts[bbin]+=1
            theta,phi = sphang(atoms[i],jatom)
            gvals[bbin]+= special.sph_harm(0,l,theta,phi)

    #At bond length 0, Qlm has one non-zero value at m=0
    Ql0 = conj(sph_harm(0,l,0,0))
    Q0=bondOrientR(atoms,basis,0,0,1) #always 0.28209479 = 1/sqrt(4*pi)
    
    #always use m=0, due to Ql0 normalizing factor which is only non-zero at m=0.
    norm  = 2*(l+1)*Q0*Q0
    for i,n in enumerate(bcnts):
        if n>0:
            w = Ql0/n/norm
            gvals[i] = (gvals[i]*w).real

    print "Finished binning bond angle values"

    return rbins,gvals


def coordinationNumber(atoms,basis,l=None,neighbs=None,rcut=None,debug=False):
    #l: not used
    
    if neighbs==None:
        if rcut==None:
            rcut = generateRCut(atoms,debug=debug)
        bounds=[[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
        neighbs = neighbors(atoms,bounds,rcut)

    cns = map(len,neighbs)
        
    return cns
            
def radialDistribution(atoms,basis,l=None,neighbs=None,rcut=None,debug=False):
     if rcut==None:
         rcut = 10.0

     rbins,rdist = rdf_periodic(atoms,basis,cutoff=rcut)
     return rbins,rdist

#translational order parameter
def translational(atoms,basis,l=None,neighbs=None,rcut=None,debug=False):
    #l: not used

    if rcut==None:
        rcut = generateRCut(atoms,debug=debug)
    r,g = rdf(atoms,cutoff=rcut)
    h=map(math.fabs,g-1)
    
    tao = sum(h)/len(r)
    return tao
        

def tetrahedral(atom,basis,l=None,neighbs=None,rcut=None,debug=False):
    pass
