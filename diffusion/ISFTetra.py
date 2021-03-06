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
from scipy.integrate import simps
from numpy import array,zeros,pi,log10,cos,sin
import pylab as pl
#mine
import utils
import parserGens,lammpsIO,outcarIO
from rootMeanSquareDist import unwrap
from vanHove import vanHoveSelf


ISFSelfSphereRefCode = """
double *qxs,*qys,*qzs,a,b,a2,b2,qx,qy,qz;
int nqs = (int)nqVecs;
qxs = (double*)malloc(sizeof(double)*(int)nqs);
qys = (double*)malloc(sizeof(double)*(int)nqs);
qzs = (double*)malloc(sizeof(double)*(int)nqs);
int qlen=0;
while( qlen < nqs ){
    a = rand()*2/(double)RAND_MAX-1.0;
    b = rand()*2/(double)RAND_MAX-1.0;
    a2 = a*a;
    b2 = b*b;

    if(a2 + b2 > 1.0) 
        continue;

    qxs[qlen] = (2*a*sqrt(1-a2-b2));
    qys[qlen] = (2*b*sqrt(1-a2-b2));
    qzs[qlen] = (1-2*(a2+b2));

    qlen++;
}

double c,d,dx,dy,dz;
double *seta,*setb;
double *setaTetra,*setbTetra;

double norm = 1.0/nAtom/nTime/nqs;
int stepSize;
for(int s=0; s<nStep; s++){ //loop over step size
  stepSize = steps[s];
  
  c=0.0;
  for(int t=0; t<nTime-stepSize; t++){ //loop over time steps
    for(int q=0;q<nqs;q++){
      qx = qxs[q];
      qy = qys[q];
      qz = qzs[q];

      seta=&(atoms[(int)(t*nAtom*3)]);
      setb=&(atoms[(int)((t+stepSize)*nAtom*3)]);
      setaTetra=&(tetra[(int)(t*nAtom)]);
      setbTetra=&(tetra[(int)((t+stepSize)*nAtom)]);

      for(int i=0; i<nAtom; i++){ //loop over atoms
        if(setbTetra[i]>0.9 || setaTetra[i]>0.9) continue;

        dx = setb[i*3+0] - seta[i*3+0];
        dy = setb[i*3+1] - seta[i*3+1];
        dz = setb[i*3+2] - seta[i*3+2];

        isfs[s] += cos(dx*qx + dy*qy + dz*qz);
   }}}
   isfs[s] *= norm;
}
"""

def ISFSelfTetra(atoms,basis,tetra,q=1.0,steps=None,nqVecs=1,kmax=1.9):#kmax (first peak=1.9, second peak=3.32)
    atoms = array(atoms)
    basis = array(basis)
    nTime = atoms.shape[0]
    nAtom = atoms.shape[1]

    if steps == None:
        nSteps = 200
        steps = sorted(list(set(array([int(10**(t*10.0/nSteps-1.0)) for t in range(nSteps)]))))
        steps.remove(0)
        steps = [t for t in steps if t<nTime]
        print steps
        steps=array(steps)

    atoms = atoms.ravel()
    tetra = tetra.ravel()

    nStep = steps.shape[0]

    isfs = zeros(nStep)
    k=kmax

    weave.inline(ISFSelfSphereRefCode,['nqVecs','steps','nStep','isfs','atoms','nTime','nAtom','k','tetra'])
    return steps,isfs

def fourierDensity(atoms,basis):
    atoms = array(atoms)
    basis = array(basis)

if __name__ == "__main__":
    utils.usage(["<dump.dat or OUTCAR>","<\'s\'-self \'d\'-distinct or \'t\'-total'>, default total"],1,3)

    RMAX = 10.0

    plotEnable = True
    logtEnable = False
    if "-noPlot" in sys.argv:
        sys.argv.remove("-noPlot")
        plotEnable = False
    if "-logt" in sys.argv:
        sys.argv.remove("-logt")
        logtEnable = True

    inputFile = sys.argv[1]
    configIterator = parserGens.parseEnsemble(inputFile+".tetra")
    tetraTime = array([tetra for tetra in configIterator])

    #steps = map(int,sys.argv[2].split(","))
    isfType = "total"
    if sys.argv[2][0] in ['s','S']:
        isfType = "self"
    else:
        isfType = "total"


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

    if isfType == "total":
        #rbins,bins = vanHoveTotal(atomsTime,basis,steps,cutr=cutr,nBin=nBin,norm=norm)
        #ylab = "$G(r,t) / V_{shell}(r)$"
        pass
    elif isfType == "self":
        steps,isfs = ISFSelfTetra(atomsTime,basis,tetraTime,nqVecs=3)
        if logtEnable:
            steps=log10(steps)
        pl.plot(steps,isfs)
        if logtEnable:
            header = "#Generated by %s using file %s\n log(t) ISF \n"%(sys.argv[0].split("/")[-1],inputFile)
            pl.xlabel("$log_{10}(t)$")
        else:
            header = "#Generated by %s using file %s\n t ISF \n"%(sys.argv[0].split("/")[-1],inputFile)
            pl.xlabel("t")
        pl.ylabel("$F_S(k,t)$")
    
    outputFile = inputFile + ".isf" + isfType[0].upper() + "tetra"
    odata = header
    for b in range(len(steps)):
        odata += str(steps[b])+" "+str(isfs[b])+"\n"
    open(outputFile,"w").write(odata)

    if plotEnable:
        pl.show()



ISFSelfGsrtRefCode = """
double c,d,dx,dy,dz;
double *seta,*setb;

int stepSize;
for(int s=0; s<nStep; s++){ //loop over step size
  stepSize = steps[s];
  
  c=0.0;
  for(int t=0; t<nTime-stepSize; t++){ //loop over time steps
    for(int q=0;q<nqs;q++){
      qx = qxs[q];
      qy = qys[q];
      qz = qzs[q];

      seta=&(atoms[(int)(t*nAtom*3)]);
      setb=&(atoms[(int)((t+stepSize)*nAtom*3)]);

      for(int i=0; i<nAtom; i++){ //loop over atoms
      
        dx = setb[i*3+0] - seta[i*3+0];
        dy = setb[i*3+1] - seta[i*3+1];
        dz = setb[i*3+2] - seta[i*3+2];

        isfs[s] += cos(dx*qx + dy*qy + dz*qz)/nAtom;
}}}}
"""

"""
def ISFSelf2(atoms,basis,q=1.0,steps=None,kmax=1.9):
    atoms = array(atoms)
    basis = array(basis)
    nTime = atoms.shape[0]
    nAtom = atoms.shape[1]

    if steps == None:
        steps = sorted(list(set(array([int(10**(t/50.0*4.0-1.0)) for t in range(50)]))))
        steps.remove(0)
        steps = [t for t in steps if t<nTime]
        steps=array(steps)

    k=kmax
    nBin=1000
    rvals,Gsrs = vanHoveSelf(atoms,basis,steps,cutr=10.0,nBin=nBin)
    for s in range(len(steps)):
        print simps([ Gsrs[s][i]*cos(rvals[i]*k) for i in range(len(rvals)) ],x=rvals)
        
    isfs = list()
    for s in range(len(steps)):
        isfs.append(simps( [ Gsrs[s][i]*cos(k*rvals[i]) for i in range(len(rvals)) ], x=rvals ))

    return steps,isfs
"""
