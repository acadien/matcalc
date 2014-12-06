#!/usr/bin/python

#Takes a dump.dat file, parses it, generates van-hove, fourier transform to
#generate intermediate scattering function

#SHOULD BE RUN ONLY ON NVT/NVE SIMULATIONS WITH UNWRAPPED COORDINATES
#Unwraps outcar coordinates

#mine
import plotRemote as pr
#theirs
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
import pylab as pl

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

      if(d/dr > nBin) continue;

      bins[s*nBin+(int)(d/dr)]+=c;
}}}
"""

def vanHoveSelf(atoms,basis,steps,cutr=10.0,nBin=2000,norm=False):
    atoms = array(atoms)
    basis = array(basis)

    #if type(steps) != type(list()): #handle the case where steps is just a single step
    #    steps = list(steps)
    #steps = array(steps)

    nTime = atoms.shape[0]
    nAtom = atoms.shape[1]
    nStep = len(steps)

    atoms = atoms.ravel()
    b = basis.ravel()

    dr = cutr/nBin
    rbins = [i*dr for i in range(nBin)]
    bins = zeros(nStep*nBin)
    weave.inline(vhSelfRefCode,['atoms','nTime','nAtom','steps','nStep','b','bins','nBin','dr'])

    bins.shape = [nStep,nBin]
    atoms.shape = [nTime,nAtom,3]
    for i in range(nStep):
        bins[i][0]=0.0

    if norm: #apply 4pir^2 correction
        for i in range(nStep):
            bins[i] = [bins[i][k]*4*pi*r*r for k,r in enumerate(rbins)] #*4*pi*r*r

    return rbins,bins

#Distinct Part of the Van Hove Function
#assumes coordinates are unwrapped
vhDistinctRefCode = """
double aix,ajx,aiy,ajy,aiz,ajz;
double c,d,dx,dy,dz;
double *seta,*setb;
double dmin;

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

        dmin=100000.;
        //Minimum image distance
        for(int t1=-1; t1<2; t1++){
        for(int t2=-1; t2<2; t2++){
        for(int t3=-1; t3<2; t3++){
            dx=aix-ajx+t1*b[0]+t2*b[3]+t3*b[6];
            dy=aiy-ajy+t1*b[1]+t2*b[4]+t3*b[7];
            dz=aiz-ajz+t1*b[2]+t2*b[5]+t3*b[8];
            d = sqrt(dx*dx + dy*dy + dz*dz);
            if(d <= dmin)
                dmin=d;
        }}}
        if( dmin < nBin*dr )
            bins[ s*nBin + (int)(dmin/dr) ]+=c;

}}}}
"""

def vanHoveDistinct(atoms,basis,steps,cutr=10.0,nBin=1000,norm=False):
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

    if norm:
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
double dmin;

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

        dmin=100000.;
        //Minimum image distance
        for(int t1=-1; t1<2; t1++){
        for(int t2=-1; t2<2; t2++){
        for(int t3=-1; t3<2; t3++){
            dx=aix-ajx+t1*b[0]+t2*b[3]+t3*b[6];
            dy=aiy-ajy+t1*b[1]+t2*b[4]+t3*b[7];
            dz=aiz-ajz+t1*b[2]+t2*b[5]+t3*b[8];
            d = sqrt(dx*dx + dy*dy + dz*dz);
            if(d <= dmin)
                dmin=d;
        }}}
        if( dmin < nBin*dr )
            bins[ s*nBin + (int)(dmin/dr) ]+=c;
}}}}
"""

def vanHoveTotal(atoms,basis,steps,cutr=10.0,nBin=1000,norm=False):
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

    if norm:
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
    utils.usage(["<dump.dat or OUTCAR>","<comma seperated step sizes [e.g. 1,2,3]>","<\'s\'-self \'d\'-distinct or \'t\'-total'>, default total"],2,4)

    RMAX = 10.0

    norm = False #an arbitrary normalization is applied
    if "-norm" in sys.argv:
        norm = True
        sys.argv.remove("-norm")

    inputFile = sys.argv[1]
    steps = map(int,sys.argv[2].split(","))
    vhType = "total"
    if len(sys.argv)==4:
        if sys.argv[3][0] in ['s','S']:
            vhType = "self"
        if sys.argv[3][0] in ['d','D']:
            vhType = "distinct"

    if "OUTCAR" in inputFile:
        ocarFile = inputFile
        atomByteNums = outcarIO.atomBytes(ocarFile)
        nAtom = outcarIO.nAtoms(ocarFile)
        basis = outcarIO.basis(ocarFile)
        
        configIterator = parserGens.parseOutcarAtoms(atomByteNums,ocarFile,nAtom)
        atomsTime = [array(atoms) for atoms in configIterator]
        atomsTime = unwrap(atomsTime,basis)
    else:
        lmpFile = inputFile
        atomByteNums = lammpsIO.atomsBytes(lmpFile)
        nAtom = lammpsIO.nAtoms(lmpFile)
        basis = lammpsIO.basis(lmpFile)
        bounds=[[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]

        configIterator = parserGens.parseLammpsAtoms(atomByteNums,lmpFile,nAtom)
        atomsTime = [array(atoms) for atoms in configIterator]

    nBin = 1000
    cutr = 10.0
    if vhType == "total":
        rbins,bins = vanHoveTotal(atomsTime,basis,steps,cutr=cutr,nBin=nBin,norm=norm)
        if norm:
            ylab = r"$G(r,t) / V_{shell}(r)$"
        else:
            ylab = r"$(G(r,t)$"
    elif vhType == "self":
        rbins,bins = vanHoveSelf(atomsTime,basis,steps,cutr=cutr,nBin=nBin,norm=norm)
        if norm:
            ylab = r"$G_s(r,t) \, $*$ 4\pi r^2 $"
        else:
            ylab = r"$G_s(r,t)$"
    elif vhType == "distinct":
        rbins,bins = vanHoveDistinct(atomsTime,basis,steps,cutr=cutr,nBin=nBin,norm=norm)
        if norm:
            ylab = r"$G_d(r,t) / V_{shell}(r)$"
        else:
            ylab = r"$G_d(r,t)$"

    xlab = "$r \; (\AA)$"
    nBin = len(rbins)

    outputFile = inputFile + ".vhove" + vhType[0].upper()
    odata = "#Generated by %s using file %s\n"%(sys.argv[0].split("/")[-1],inputFile)
    odata += r"$r(\AA)$ " + "ts ".join(map(str,steps))+"ts\n"
    transposeBins=zip(*bins)
    for b in range(nBin):
        odata += str(rbins[b])+" "+" ".join(map(str,transposeBins[b]))+"\n"
    open(outputFile,"w").write(odata)

    for i,s in enumerate(steps):
        pl.plot(rbins,bins[i],label=str(steps[i]))
    pl.legend(loc=0,title="step size")
    pl.xlabel(xlab)
    pl.ylabel(ylab)
    pl.title("VanHove - "+vhType)
    pr.prshow()


