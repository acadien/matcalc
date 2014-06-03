#!/usr/bin/python
#Raw functions for calculating order parameters such as bond orientation, translational and tetrahedral ordering.

import numpy as np
from scipy.special import sph_harm as sph_harm
from scipy import conj,array,integrate
import scipy
import math
from operator import mul
import pylab as pl
import plotRemote as pr
#mine
from neighbors import neighbors,nNearestNeighbors,full2half,voronoiNeighbors
from struct_tools import *
from rdf import *
from sf import sf

#local bond-orientational: Q_l for each atom.  atomi>-1 selects a specific atom
#rcut is interms of shells not a distance
def bondOrientation(atoms,basis,l,neighbs=None,rcut=1,debug=False):

    if neighbs==None:
        bounds=[[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
        if rcut<=1:
            rcut = generateRCut(atoms,basis,debug=debug)
            print "Automatically generating r-cutoff=",rcut
            neighbs = neighbors(atoms,bounds,rcut)
        elif rcut==2:
            rcut = generateRCut(atoms,basis,debug=debug)
            print "Automatically generating r-cutoff=",rcut
            neighbs = neighbors(atoms,bounds,rcut)
            neighbs = secondShell(neighbs)
        else:
            neighbs = neighbors(atoms,bounds,rcut)

    #sum the spherical harmonic over ever neighbor pair
    Qlms = [sum( [ pairSphereHarms(atoms[i],minImageAtom(atoms[i],atoms[j],basis),l) for j in ineighbs ] ) / len(ineighbs) for i,ineighbs in enumerate(neighbs) ] 
    Ql = [ (((Qlm.conjugate()*Qlm *4*np.pi / (2*l+1.))).real)**0.5 for Qlm in Qlms] 

    return Ql,rcut

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
    if rcut==None:
        rcut = generateRCut(atoms,basis,debug=debug)

    if neighbs==None:
        rcut = 6.0
        bounds=[[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
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
            rcut = generateRCut(atoms,basis,debug=debug)
            print "Using RDF to generate r-cutoff=",rcut
        else:
            "Using r-cutoff=",rcut

        bounds=[[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
        neighbs = neighbors(atoms,bounds,rcut,style="full")
        #neighbs = voronoiNeighbors(atoms,basis,[1]*len(atoms),style="full")
    cns = map(len,neighbs)
        
    return cns,rcut

def radangDistribution(atoms,basis,l=None,neighbs=None,rcut=None,debug=False):
    #l: not used
    if neighbs==None:
        if rcut==None:
            rcut = generateRCut(atoms,basis,debug=debug)
            print "Using RDF to generate r-cutoff=",rcut
        else:
            print "Using r-cutoff=",rcut

        bounds=[[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
        neighbs = neighbors(atoms,bounds,rcut,style="full")
    return rdf_by_adf(atoms,neighbs,basis,rcut=rcut)
            
def radialDistribution(atoms,basis,l=1000,neighbs=None,rcut=None,debug=False):
     if rcut==None:
         rcut = 10.0
 
     nbins=l
     if nbins==None:
         nbins=1000

     return rdf_periodic(atoms,basis,cutoff=rcut,nbins=nbins)#,rdist

def angleDistribution(atoms,basis,l=None,neighbs=None,rcut=None,debug=False):
    if rcut==None:
        rcut = generateRCut(atoms,basis,debug=debug)
        print "Using RDF to generate r-cutoff=",rcut
    else:
        print "Using r-cutoff=",rcut

    if neighbs==None:
        bounds = [[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
        #neighbs = voronoiNeighbors(atoms,basis,[1]*len(atoms),style="full")
        neighbs = neighbors(atoms,bounds,rcut)

    return adf(atoms,neighbs,basis,rcut,nbins=360)

def structureFactor(atoms,basis,l=None,neighbs=None,rcut=None,debug=False):
    if rcut==None:
        rcut = min(sum([basis[0][0],basis[1][1],basis[2][2]])/6, 10.0)
        print "Automatically generating r-cutoff=",rcut

    if l==None:
        l=12.0

    rbins,rdist = rdf_periodic(atoms,basis,cutoff=rcut)
    Nr=len(rbins)
    density = atoms.shape[0] / volume(basis)

    qbins,qvals = sf(rbins,rdist,density,Lmax=l,qbins=2048,damped=True)
    return qbins,qvals

#def structureFactor0(atoms,basis,l=None,neighbs=None,rcut=None,debug=False):
#   sfq(atoms,basis)
#   return 0
            
def structureFactor0(rdfX,rdfY,nDensity,l=None,neighbs=None,rcut=None,debug=False):
    rdfX=array(rdfX)
    rdfY=array(rdfY)
    hr = rdfY-1.0
    from scipy import fftpack
    sp = fftpack.rfft(hr*rdfX*rdfX)
    print sp
    freq = fftpack.fftfreq(len(hr),d=(rdfX[1]-rdfX[0]))

    import pylab as pl
    pl.plot(freq,((sp.real/288.)**2+(sp.imag/288)**2)**0.5+1.0)
    pl.plot(freq,sp.real/288.+1.0,freq,sp.imag/288+1.0)
    pl.show()
    exit(0)
    return [integrate.simps((rdfY-1.0),rdfX)*nDensity+1]
    
#translational order parameter, l=neighbor shell
def translational(atoms,basis,l=None,neighbs=None,rcut=None,debug=False):
    #l: not used

    if rcut==None:
        rcut = 10.0#generateRCut(atoms,basis,debug=debug)
    r,g = rdf_periodic(atoms,basis,cutoff=rcut)#rbins,rdist

    h=map(math.fabs,g-1)

    tao = [i/len(r) for i in h]
    return tao,rcut
        
#Sg and Sk as defined by:
#P.L. Chau and A.J. Hardwick, J. Mol. Phys, V93, pp511-518, No3, (1998)
#Sg is 1 for tetrahedral, 3/4 for randomly arranged bonds and 0 for superimposed
def tetrahedral(atoms,basis,l=None,neighbs=None,rcut=None,debug=False):
    atoms=array(atoms)

    if neighbs==None:
        if rcut==None:
            rcut = generateRCut(atoms,basis,debug=debug)
            rcut += 2.5
        bounds=[[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
        #ensure only the 4 shortest bonds are used
        neighbs = nNearestNeighbors(4,atoms,bounds,rcut)
            
    third=1./3.
    tets=list()
    for i,ineighbs in enumerate(neighbs):
        iatom=atoms[i]

        if len(ineighbs)<3:
            tets.append(0)
            continue

        Sg=0
        for v,j in enumerate(ineighbs):
            jatom=minImageAtom(iatom,atoms[j],basis)

            for k in ineighbs[v+1:]:
                katom=minImageAtom(iatom,atoms[k],basis)
                a = ang(iatom,jatom,katom)
                Sg+=(cos(a)+third)**2
        Sg = 1 - 3./8.*Sg
        tets.append(Sg)
    return tets,rcut
