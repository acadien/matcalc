#!/usr/bin/python

#Takes a dump.dat file, parses it, generates van-hove, fourier transform to
#generate intermediate scattering function

#SHOULD BE RUN ONLY ON NVT/NVE SIMULATIONS WITH UNWRAPPED COORDINATES
#Unwraps outcar coordinates

import sys
import itertools
from scipy import weave
from scipy.weave import converters
from struct_tools import dist
from numpy import array,zeros,pi
#mine
import utils
import parserGens,lammpsIO,outcarIO
from rootMeanSquareDist import unwrap
from struct_tools import volume

#Self Part of the Van Hove Function
vhSelfRefCode = """
double c,d,dx,dy,dz;
double *seta,*setb;

int stepSize;
for(int s=0; s<nStep; s++){ //loop over step size
  stepSize = steps[s];
  
  c = 1.0/(nTime-stepSize)/nAtom;
  for(int t=0; t<nTime-stepSize; t++){ //loop over time steps
    seta=&(atoms[(int)(t*nAtom*3)]);
    setb=&(atoms[(int)((t+stepSize)*nAtom*3)]);

    for(int i=0; i<nAtom; i++){ //loop over atoms
      
      dx = setb[i*3+0] - seta[i*3+0];
      dy = setb[i*3+1] - seta[i*3+1];
      dz = setb[i*3+2] - seta[i*3+2];
      d = sqrt(dx*dx + dy*dy + dz*dz);

      bins[s*nBin+(int)(d/dr)]+=c;
}}}
"""

def vanHoveSelf(atoms,basis,steps,cutr=10.0,nBin=1000):
    atoms = array(atoms)
    basis = array(basis)
    
    if type(steps) != type(list()): #handle the case where steps is just a single step
        steps = list(steps)
    steps = array(steps)

    nTime = atoms.shape[0]
    nAtom = atoms.shape[1]
    nStep = steps.shape[0]

    atoms = atoms.ravel()
    b = basis.ravel()

    dr = cutr/nBin
    rbins = [i*dr for i in range(nBin)]
    bins = zeros(nStep*nBin)
    weave.inline(vhSelfRefCode,['atoms','nTime','nAtom','steps','nStep','b','bins','nBin','dr'])

    bins.shape = [nStep,nBin]
    atoms.shape = [nTime,nAtom,3]

    return rbins,bins

#Distinct Part of the Van Hove Function
#assumes coordinates are unwrapped
vhDistinctRefCode = """
double aix,ajx,aiy,ajy,aiz,ajz;
double c,d,dx,dy,dz;
double *seta,*setb;
double cmin;

int stepSize;
for(int s=0; s<nStep; s++){ //loop over step size
  stepSize = steps[s];
  
  c = 2.0/(nTime-stepSize)/nAtom;
  for(int t=0; t<nTime-stepSize; t++){ //loop over time steps
    seta=&(atoms[(int)(t*nAtom*3)]);
    setb=&(atoms[(int)((t+stepSize)*nAtom*3)]);

    d=0;
    for(int i=0; i<nAtom; i++){ //loop over atoms
      aix = seta[i*3];
      aiy = seta[i*3+1];
      aiz = seta[i*3+2];
     
      for(int j=i+1; j<nAtom; j++){ //loop over other atoms
        ajx = setb[j*3];
        ajy = setb[j*3+1];
        ajz = setb[j*3+2];

        cmin=100000.;
        //Minimum image distance
        for(int t1=-1; t1<2; t1++){
        for(int t2=-1; t2<2; t2++){
        for(int t3=-1; t3<2; t3++){
            dx=aix-ajx+t1*b[0]+t2*b[3]+t3*b[6];
            dy=aiy-ajy+t1*b[1]+t2*b[4]+t3*b[7];
            dz=aiz-ajz+t1*b[2]+t2*b[5]+t3*b[8];
            d = sqrt(dx*dx + dy*dy + dz*dz);
            if(d <= cmin)
                cmin=d;
        }}}
        if( cmin < nBin*dr )
            bins[ s*nBin + (int)(cmin/dr) ]+=c;

}}}}

"""

def vanHoveDistinct(atoms,basis,steps,cutr=10.0,nBin=1000):
    atoms = array(atoms)
    basis = array(basis)
    
    if type(steps) != type(list()): #handle the case where steps is just a single step
        steps = list(steps)
    steps = array(steps)

    nTime = atoms.shape[0]
    nAtom = atoms.shape[1]
    nStep = steps.shape[0]

    atoms = atoms.ravel()
    b = basis.ravel()

    dr = cutr/nBin
    rbins = [i*dr for i in range(nBin)]
    bins = zeros(nStep*nBin)
    weave.inline(vhDistinctRefCode,['atoms','nTime','nAtom','steps','nStep','b','bins','nBin','dr'])

    bins.shape = [nStep,nBin]
    atoms.shape = [nTime,nAtom,3]

    Ndensity=nAtom/volume(basis)
    for i,r in enumerate(rbins):
        if i==0:
            vol=4.0*pi*dr*dr*dr/3.0
        else:
            vol=4.0*pi*r*r*dr
        for k in range(nStep):
            bins[k][i] /= vol*Ndensity

    return rbins,bins


#Total Van Hove function (sum of self and distinct)
#assumes coordinates are unwrapped
vhTotalRefCode = """
double aix,ajx,aiy,ajy,aiz,ajz;
double c,d,dx,dy,dz;
double *seta,*setb;
double cmin;

int stepSize;
for(int s=0; s<nStep; s++){ //loop over step size
  stepSize = steps[s];
  
  c = 2.0/(nTime-stepSize)/nAtom;
  for(int t=0; t<nTime-stepSize; t++){ //loop over time steps
    seta=&(atoms[(int)(t*nAtom*3)]);
    setb=&(atoms[(int)((t+stepSize)*nAtom*3)]);

    d=0;
    for(int i=0; i<nAtom; i++){ //loop over atoms
      aix = seta[i*3];
      aiy = seta[i*3+1];
      aiz = seta[i*3+2];
     
      for(int j=i; j<nAtom; j++){ //loop over other atoms
        ajx = setb[j*3];
        ajy = setb[j*3+1];
        ajz = setb[j*3+2];

        cmin=100000.;
        //Minimum image distance
        for(int t1=-1; t1<2; t1++){
        for(int t2=-1; t2<2; t2++){
        for(int t3=-1; t3<2; t3++){
            dx=aix-ajx+t1*b[0]+t2*b[3]+t3*b[6];
            dy=aiy-ajy+t1*b[1]+t2*b[4]+t3*b[7];
            dz=aiz-ajz+t1*b[2]+t2*b[5]+t3*b[8];
            d = sqrt(dx*dx + dy*dy + dz*dz);
            if(d <= cmin)
                cmin=d;
        }}}
        if( cmin < nBin*dr )
            bins[ s*nBin + (int)(cmin/dr) ]+=c;
}}}}
"""

def vanHoveTotal(atoms,basis,steps,cutr=10.0,nBin=1000):
    atoms = array(atoms)
    basis = array(basis)
    
    if type(steps) != type(list()): #handle the case where steps is just a single step
        steps = list(steps)
    steps = array(steps)

    nTime = atoms.shape[0]
    nAtom = atoms.shape[1]
    nStep = steps.shape[0]

    atoms = atoms.ravel()
    b = basis.ravel()

    dr = cutr/nBin
    rbins = [i*dr for i in range(nBin)]
    bins = zeros(nStep*nBin)
    weave.inline(vhTotalRefCode,['atoms','nTime','nAtom','steps','nStep','b','bins','nBin','dr'])

    bins.shape = [nStep,nBin]
    atoms.shape = [nTime,nAtom,3]

    Ndensity=nAtom/volume(basis)
    for i,r in enumerate(rbins):
        if i==0:
            vol=4.0*pi*dr*dr*dr/3.0
        else:
            vol=4.0*pi*r*r*dr
        for k in range(nStep):
            bins[k][i] /= vol*Ndensity

    return rbins,bins


if __name__ == "__main__":
    utils.usage(["<dump.dat or OUTCAR>","<comma seperated step sizes [e.g. 1,2,3]>","<\'s\'-self \'d\'-distinct or \'t\'-total'>, default total"],2,3)

    RMAX = 10.0

    inputFile = sys.argv[1]
    steps = map(int,sys.argv[2].split(","))
    vhType = "total"
    if len(sys.argv)==4:
        if sys.argv[3][0] in ['s','S']:
            vhType = "self"
        if sys.argv[3][0] in ['d','D']:
            vhType = "distinct"

    if "OUTCAR" in inputFile:
        pass
        unwrap(atoms,basis)
    else:
        lmpFile = inputFile
        atomByteNums = lammpsIO.atomsBytes(lmpFile)
        nAtom = lammpsIO.nAtoms(lmpFile)
        basis = lammpsIO.basis(lmpFile)
        bounds=[[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]

        nStep=len(steps)

        configIterator = parserGens.parseLammpsAtoms(atomByteNums,lmpFile,nAtom)
        atomsTime = [array(atoms) for atoms in configIterator]


    nBin = 1000
    cutr = 10.0

    if vhType == "total":
        rbins,bins = vanHoveTotal(atomsTime,basis,steps,cutr=cutr,nBin=nBin)
        ylab = "$G(r,t) / V_{shell}(r)$"
    elif vhType == "self":
        rbins,bins = vanHoveSelf(atomsTime,basis,steps,cutr=cutr,nBin=nBin)
        for i in range(nStep):
            bins[i] = [bins[i][k]*4*pi*rbins[k]**2 for k in range(nBin)]
        ylab = "$4\pi r^2 $*$ \, G_s(r,t)$"
    elif vhType == "distinct":
        rbins,bins = vanHoveDistinct(atomsTime,basis,steps,cutr=cutr,nBin=nBin)
        ylab = r"$G_d(r,t) / V_{shell}(r)$"

    xlab = "$r \; (\AA)$"
    nBin = len(rbins)

    import pylab as pl
    for i,s in enumerate(steps):
        pl.plot(rbins,bins[i],label=str(steps[i]))
    pl.legend(loc=0,title="step size")
    pl.xlabel(xlab)
    pl.ylabel(ylab)
    pl.title("VanHove - "+vhType)
    pl.show()
