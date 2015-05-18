#!/usr/bin/python
#Raw functions for calculating order parameters such as bond orientation, translational and tetrahedral ordering.

import numpy as np
from scipy.special import sph_harm as sph_harm
from scipy import conj,array,integrate,weave
from scipy.weave import converters
import scipy
import math
from operator import mul
import pylab as pl
import plotRemote as pr
#mine
from neighbors import neighbors,nNearestNeighbors,full2half,voronoiNeighbors,secondShell
from struct_tools import *
from rdf import *
from sf import sf,sfq,sfq0,sfq2

#Maps atomic coordinates back into the basis, assumes an orthogonal basis set
def rectify(atoms,basis):
    atoms[:,0] = np.mod(atoms[:,0],basis[0][0])
    atoms[:,1] = np.mod(atoms[:,1],basis[1][1])
    atoms[:,2] = np.mod(atoms[:,2],basis[2][2])
    return atoms

#local bond-orientational: Q_l for each atom.  atomi>-1 selects a specific atom
#rcut is interms of shells not a distance
def bondOrientation(atoms,basis,l,neighbs=None,rcut=None,debug=False):
    atoms = array(atoms)
    basis = array(basis)    
    atoms = rectify(atoms,basis)

    if neighbs==None:
        bounds=[[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
        
        if rcut==None:
            rcut = generateRCut(atoms,basis,debug=debug)
            #print "Automatically generating r-cutoff=",rcut
        
        neighbs = neighbors(atoms,bounds,rcut)

    #sum the spherical harmonic over ever neighbor pair
    a = 4*np.pi / (2*l+1.)
    Ql=list()
    for i,ineighbs in enumerate(neighbs):
        n=len(ineighbs)

        shij = np.vectorize(complex)(zeros(2*l+1)) #spherical harmonic for bond i-j
        for j in ineighbs:
            shij += pairSphereHarms(atoms[i],minImageAtom(atoms[i],atoms[j],basis),l)/n
        shi = a * sum( scipy.real( scipy.multiply(shij,scipy.conj(shij)) ) )
        Ql.append(shi**0.5)
    return Ql,rcut

#local bond-orientational: Q_l for each atom.  atomi>-1 selects a specific atom
#rcut is interms of shells not a distance
#uses 2nd shell from neighbor list
def bondOrientation2sh(atoms,basis,l,neighbs=None,rcut=None,debug=False):
    atoms = array(atoms)
    basis = array(basis)    
    atoms = rectify(atoms,basis)

    if neighbs==None:
        bounds=[[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]

        if rcut==None:
            rcut = generateRCut(atoms,basis,debug=debug)
            #print "Automatically generating r-cutoff=",rcut

        neighbs = secondShell( neighbors(atoms,bounds,rcut) )

    #sum the spherical harmonic over ever neighbor pair
    a = 4*np.pi / (2*l+1.)
    Ql=list()
    for i,ineighbs in enumerate(neighbs):
        n=len(ineighbs)

        shij = np.vectorize(complex)(zeros(2*l+1)) #spherical harmonic for bond i-j
        for j in ineighbs:
            shij += pairSphereHarms(atoms[i],minImageAtom(atoms[i],atoms[j],basis),l)/n
        shi = a * sum( scipy.real( scipy.multiply(shij,scipy.conj(shij)) ) )
        Ql.append(shi**0.5)
    
    return Ql,rcut

#Helper function, returns Qlm values m=(-l .. 0 .. +l) for a specific atom pair: atomi,atomj
def bondOrientR(atoms,basis,l,atomi,atomj):
    atoms = array(atoms)
    basis = array(basis)    
    atoms = rectify(atoms,basis)

    ai = atoms[atomi]
    aj = atoms[atomj]
    Qlm = pairSphereHarms(ai,minImageAtom(ai,aj,basis),l)
    return Qlm

#Bond antle correlation function Gl as defined in:
#       Nature Materials, Vol2, Nov. 2003, Sastry & Angell
def bondAngleCorr(atoms,basis,l,neighbs=None,rcut=None,debug=False):
    atoms = array(atoms)
    basis = array(basis)    
    atoms = rectify(atoms,basis)

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
    atoms = array(atoms)
    basis = array(basis)    
    atoms = rectify(atoms,basis)

    #l: not used
    if neighbs==None:
        if rcut==None:
            rcut = generateRCut(atoms,basis,debug=debug)
            #print "Using RDF to generate r-cutoff=",rcut
        else:
            "Using r-cutoff=",rcut

        bounds=[[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
        neighbs = neighbors(atoms,bounds,rcut,style="full")
        #neighbs = voronoiNeighbors(atoms,basis,[1]*len(atoms),style="full")
    cns = map(len,neighbs)
        
    return cns,rcut

def abraham(atoms,basis,l=None,neighbs=None,rcut=None,debug=False):
    atoms = array(atoms)
    basis = array(basis)    
    atoms = rectify(atoms,basis)

    #Generate g(r)
    rvals,gr = rdf_periodic(atoms,basis,cutoff=6.0)
    
    #Smoothed g(r)
    sgr=windowAvg(gr,n=25).tolist()
    #derivative of smoothed-g(r)
    dsgr = windowAvg(windowAvg([(sgr[i+1]-sgr[i])/(rvals[1]-rvals[0]) for i in range(len(sgr)-1)],n=50),n=20).tolist()

    first_neg = [i for i,v in enumerate(dsgr) if v<0][0]
    first_peak = sgr.index(max(sgr[:first_neg]))
    rindex = first_neg + [i for i,v in enumerate(dsgr[first_neg:]) if v>=0][0]
    rmin = rvals[rindex]
    gmin = sgr[rindex]
    rmax = rvals[first_peak-2]
    gmax = sgr[first_peak-2]

    return gmin/gmax

def radangDistribution(atoms,basis,l=None,neighbs=None,rcut=None,debug=False):
    atoms = array(atoms)
    basis = array(basis)    
    atoms = rectify(atoms,basis)

    #l: not used
    if neighbs==None:
        if rcut==None:
            rcut = generateRCut(atoms,basis,debug=debug)
            #print "Using RDF to generate r-cutoff=",rcut
        #else:
        #    print "Using r-cutoff=",rcut

        bounds=[[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
        neighbs = neighbors(atoms,bounds,rcut,style="full")
    return rdf_by_adf(atoms,neighbs,basis,rcut=rcut)
            
def radialDistribution(atoms,basis,l=1000,neighbs=None,rcut=None,debug=False):
    atoms = array(atoms)
    basis = array(basis)    
    atoms = rectify(atoms,basis)

    if rcut==None:
        rcut = 10.0
        
    nbins=l
    if nbins==None:
        nbins=1000
     
    return rdf_periodic(atoms,basis,cutoff=rcut,nbins=nbins)#,rdist

def angleDistribution(atoms,basis,l=None,neighbs=None,rcut=None,debug=False):
    atoms = array(atoms)
    basis = array(basis)    
    atoms = rectify(atoms,basis)

    if rcut==None:
        rcut = generateRCut(atoms,basis,debug=debug)
    #    print "Using RDF to generate r-cutoff=",rcut
    #else:
    #    print "Using r-cutoff=",rcut

    if neighbs==None:
        bounds = [[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
        #neighbs = voronoiNeighbors(atoms,basis,[1]*len(atoms),style="full")
        neighbs = neighbors(atoms,bounds,rcut)

    return adf(atoms,neighbs,basis,rcut,nbins=360)

def structureFactor(atoms,basis,l=None,neighbs=None,rcut=None,debug=False):
    atoms = array(atoms)
    basis = array(basis)    
    atoms = rectify(atoms,basis)

    #if rcut==None:
    #    rcut = min(sum([basis[0][0],basis[1][1],basis[2][2]])/6, 10.0)
    #    print "Automatically generating r-cutoff=",rcut

    if l==None:
        l=12.0

    #rbins,rdist = rdf_periodic(atoms,basis,cutoff=rcut)
    #Nr=len(rbins)
    #density = atoms.shape[0] / volume(basis)

    #qbins,qvals = sf(rbins,rdist,density,Lmax=l,qbins=2048,damped=False)
    qbins,qvals = sfq(atoms,basis,nqbins=290,qcut=l)
    return qbins,qvals

def structureFactor0(atoms,basis,l=None,neighbs=None,rcut=None,debug=False):
    atoms = array(atoms)
    basis = array(basis)    
    atoms = rectify(atoms,basis)

    #sfq(atoms,basis)
    if rcut==None:
        rcut = sum([basis[0][0],basis[1][1],basis[2][2]])/6
        if rcut<20.0: rcut = 10.0
    #    print "Automatically generating r-cutoff=",rcut

    if l==None:
        l=6.0
    
    rcut = float(rcut)
    rbins,rdist = rdf_periodic(atoms,basis,cutoff=rcut)
    Nr=len(rbins)
    density = atoms.shape[0] / volume(basis)

    qbins,qvals = sfq0(rbins,rdist,density,Lmax=l,qbins=2048,damped=True)

    return 0

#translational order parameter, l=neighbor shell
def translational(atoms,basis,l=None,neighbs=None,rcut=None,debug=False):
    #l: not used

    if rcut==None:
        rcut = 10.0#generateRCut(atoms,basis,debug=debug)
    r,g = rdf_periodic(atoms,basis,cutoff=rcut)#rbins,rdist

    h=map(math.fabs,g-1)

    tao = [i/len(r) for i in h]
    return tao,rcut


tetraCode="""
//For mirroring atoms
double aix,aiy,aiz,ajx,ajy,ajz,akx,aky,akz;
double dx,dy,dz,c,minc;
int mint1,mint2,mint3,j,k,jzero;

//For angle calculation
double acosAng,dkx,dky,dkz,djx,djy,djz,djr,dkr;

//For Sg calculation
double third=1.0/3.0,Sg;

for(int i=0;i<nAtoms;i++){ //Central atom
  aix = atoms[i*3+0];
  aiy = atoms[i*3+1];
  aiz = atoms[i*3+2];

  Sg=0.0;

  if(i==0) jzero = 0;
  else jzero = nNeighbs[i-1];

  for(int jj=jzero;jj<nNeighbs[i];jj++){
    j = neighbs[jj];
    ajx = atoms[j*3+0];
    ajy = atoms[j*3+1];
    ajz = atoms[j*3+2];

    //Periodic Bounds - Minimum image atom between j and i
    minc = 10000.0;
    for(int t1=-1;t1<2;t1++){
    for(int t2=-1;t2<2;t2++){
    for(int t3=-1;t3<2;t3++){
      dx = ajx + t1*b[0]+t2*b[3]+t3*b[6] - aix;
      dy = ajy + t1*b[1]+t2*b[4]+t3*b[7] - aiy;
      dz = ajz + t1*b[2]+t2*b[5]+t3*b[8] - aiz;

      c = dx*dx+dy*dy+dz*dz;

      if( c <= minc ){
        mint1 = t1; mint2 = t2; mint3 = t3;
        minc = c;
      }
    }}}
    ajx += mint1*b[0] + mint2*b[3] + mint3*b[6];
    ajy += mint1*b[1] + mint2*b[4] + mint3*b[7];
    ajz += mint1*b[2] + mint2*b[5] + mint3*b[8];
    
  for(int kk=jj+1;kk<nNeighbs[i];kk++){
    k = neighbs[kk];
    akx = atoms[k*3+0];
    aky = atoms[k*3+1];
    akz = atoms[k*3+2];

    //Periodic Bounds - Minimum image atom between k and i
    minc = 10000.0;
    for(int t1=-1;t1<2;t1++){
    for(int t2=-1;t2<2;t2++){
    for(int t3=-1;t3<2;t3++){
      dx = akx + t1*b[0]+t2*b[3]+t3*b[6] - aix;
      dy = aky + t1*b[1]+t2*b[4]+t3*b[7] - aiy;
      dz = akz + t1*b[2]+t2*b[5]+t3*b[8] - aiz;

      c = dx*dx+dy*dy+dz*dz;

      if( c <= minc ){
        mint1 = t1; mint2 = t2; mint3 = t3;
        minc = c;
      }
    }}}
    akx += mint1*b[0] + mint2*b[3] + mint3*b[6];
    aky += mint1*b[1] + mint2*b[4] + mint3*b[7];
    akz += mint1*b[2] + mint2*b[5] + mint3*b[8];

    //Calculate the angle between i,j,k
    dkx = akx - aix;
    dky = aky - aiy;
    dkz = akz - aiz;
    djx = ajx - aix;
    djy = ajy - aiy;
    djz = ajz - aiz;
    djr = djx*djx + djy*djy + djz*djz;
    dkr = dkx*dkx + dky*dky + dkz*dkz;
    acosAng = (djx*dkx + djy*dky + djz*dkz)/sqrt(djr*dkr);

    c = acosAng + third;

    Sg += c*c;

  }}

  tets[i] = 1 - 3.*Sg/8.;

}
"""


#Sg and Sk as defined by:
#P.L. Chau and A.J. Hardwick, J. Mol. Phys, V93, pp511-518, No3, (1998)
#Sg is 1 for tetrahedral, 3/4 for randomly arranged bonds and 0 for superimposed
def tetrahedral(atoms,basis,l=None,neighbs=None,rcut=None,debug=False):
    atoms = array(atoms)
    basis = array(basis)
    atoms = rectify(atoms,basis)

    if neighbs==None:
        if rcut==None:
            rcut = generateRCut(atoms,basis,debug=debug)
            rcut += 3.0
        bounds=[[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
        #ensure only the 4 shortest bonds are used
        neighbs = nNearestNeighbors(4,atoms,bounds,rcut)

    def accumulateLen(x, y=[0]): y[0] += len(x); return y[0];

    nNeighbs = array(map(accumulateLen,neighbs))
    nAtoms = atoms.shape[0]
    atoms.shape = nAtoms*3
    neighbs = concatenate(map(array,neighbs))

    tets=zeros(nAtoms)    
    b=basis
    b.shape=9
    weave.inline(tetraCode,['tets','atoms','nAtoms','neighbs','nNeighbs','b'])
    atoms.shape=[nAtoms,3]

    return tets,rcut
