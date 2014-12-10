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
from numpy import *
import pylab as pl
#mine
import utils
import parserGens,lammpsIO,outcarIO
from rootMeanSquareDist import unwrap
from vanHove import vanHoveSelf

ISFFullSphereRefCode="""
//Calculate Fourier Density at each time step
double *qxs,*qys,*qzs,a,b,a2,b2,qx,qy,qz;
double *rhoqcos,*rhoqsin,*seta,*setb;
double da,db;
int stepSize;
qxs = (double*)malloc(sizeof(double)*(int)nqVecs);
qys = (double*)malloc(sizeof(double)*(int)nqVecs);
qzs = (double*)malloc(sizeof(double)*(int)nqVecs);
rhoqcos = (double*)calloc(nTime,sizeof(double));
rhoqsin = (double*)calloc(nTime,sizeof(double));

int qlen=0;
while( qlen < nqVecs ){
  a = rand()*2/(double)RAND_MAX-1.0;
  b = rand()*2/(double)RAND_MAX-1.0;
  a2 = a*a;
  b2 = b*b;

  if(a2 + b2 > 1.0) 
    continue;
  qlen++;

  qx = (2*a*sqrt(1-a2-b2))*q;
  qy = (2*b*sqrt(1-a2-b2))*q;
  qz = (1-2*(a2+b2))*q;

  for(int s=0; s<nStep; s++){
    stepSize = steps[s];
    for(int t=0; t<nTime-stepSize; t++){
      seta = &(atoms[(int)(t*nAtom*3)]);
      setb = &(atoms[(int)((t+stepSize)*nAtom*3)]);
      
      for(int i=0; i<nAtom; i++){
        da = seta[i*3+0]*qx + seta[i*3+1]*qy + seta[i*3+2]*qz;
        db = setb[i*3+0]*qx + setb[i*3+1]*qy + setb[i*3+2]*qz;
      
        isfs[s] += cos(da)*cos(db) + sin(da)*sin(db);
      }
    }
    isfs[s] /= (nTime-stepSize)*nAtom;
  }
}
"""

def ISFFullRmsd(atoms,basis,rmsd,steps=None,nqVecs=3,qmax=3.32,nStep=250,criteria=None):#qmax (first peaq=1.9, second peaq=3.32)
    atoms = array(atoms)
    basis = array(basis)
    nTime = atoms.shape[0]
    nAtom = atoms.shape[1]

    steps = sorted(list(set(array([int(10**(t*10.0/nStep-1.0)) for t in range(nStep)]))))
    steps.remove(0)
    steps = [t for t in steps if t<nTime]
    steps = array(steps)
    nStep = steps.shape[0]

    nDel=0
    ltFlag = False
    if "lt" in criteria:
        ltFlag=True
    cutoff = float(criteria[2:])

    for i,r in enumerate(rmsd):
        if ltFlag:
            if r>cutoff:
                atoms = delete(atoms,i-nDel,1)
                nAtom -= 1
                nDel +=1
        else:
            if r<cutoff:
                atoms = delete(atoms,i-nDel,1)
                nAtom -= 1
                nDel +=1
    atoms = atoms.ravel()
    isfs = zeros(nStep)
    q=qmax
    
    weave.inline(ISFFullSphereRefCode,['nqVecs','steps','nStep','isfs','atoms','nTime','nAtom','q'])

    return steps,isfs

ISFSelfSphereRefCode = """
double *qxs,*qys,*qzs,a,b,a2,b2,qx,qy,qz;
qxs = (double*)malloc(sizeof(double)*(int)nqVecs);
qys = (double*)malloc(sizeof(double)*(int)nqVecs);
qzs = (double*)malloc(sizeof(double)*(int)nqVecs);
int qlen=0;
while( qlen < nqVecs ){
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
double norm = 1.0/nAtom/nTime/nqVecs;
int stepSize;
for(int s=0; s<nStep; s++){ //loop over step size
  stepSize = steps[s];
  
  c=0.0;
  for(int t=0; t<nTime-stepSize; t++){ //loop over time steps
    for(int q=0;q<nqVecs;q++){
      qx = qxs[q];
      qy = qys[q];
      qz = qzs[q];

      seta=&(atoms[(int)(t*nAtom*3)]);
      setb=&(atoms[(int)((t+stepSize)*nAtom*3)]);

      for(int i=0; i<nAtom; i++){ //loop over atoms

        dx = setb[i*3+0] - seta[i*3+0];
        dy = setb[i*3+1] - seta[i*3+1];
        dz = setb[i*3+2] - seta[i*3+2];

        isfs[s] += cos(dx*qx + dy*qy + dz*qz);
   }}}
   isfs[s] *= norm;
}
"""

def ISFSelfRmsd(atoms,basis,rmsd,nqVecs=1,qmax=3.32,nStep=250,criteria=None):#qmax (first peaq=1.9, second peaq=3.32)
    atoms = array(atoms)
    basis = array(basis)
    nTime = atoms.shape[0]
    nAtom = atoms.shape[1]

    steps = sorted(list(set(array([int(10**(t*10.0/nStep-1.0)) for t in range(nStep)]))))
    steps.remove(0)
    steps = [t for t in steps if t<nTime]
    steps = array(steps)
    nStep = len(steps)

    nDel=0
    ltFlag = False
    if "lt" in criteria:
        ltFlag=True
    cutoff = float(criteria[2:])

    for i,r in enumerate(rmsd):
        if ltFlag:
            if r>cutoff:
                atoms = delete(atoms,i-nDel,1)
                nAtom -= 1
                nDel +=1
        else:
            if r<cutoff:
                atoms = delete(atoms,i-nDel,1)
                nAtom -= 1
                nDel +=1

    atoms = atoms.ravel()

    nStep = steps.shape[0]

    isfs = zeros(nStep)
    q=qmax

    weave.inline(ISFSelfSphereRefCode,['nqVecs','steps','nStep','isfs','atoms','nTime','nAtom','q','rmsd'])

    return steps,isfs


if __name__ == "__main__":
    utils.usage(["<dump.dat or OUTCAR>","<\'s\'-self or \'t\'-total'>,  <critera e.g. lt0.8, gt1.3>"],1,8)

    RMAX = 10.0

    plotEnable = True
    logtEnable = False
    scale = None #units in seconds
    nStep = 250
    if "-noPlot" in sys.argv:
        sys.argv.remove("-noPlot")
        plotEnable = False
    if "-logt" in sys.argv:
        sys.argv.remove("-logt")
        logtEnable = True
    if "-scale" in sys.argv:
        i = sys.argv.index("-scale")
        scale = float(sys.argv[i+1])
        sys.argv.pop(i)
        sys.argv.pop(i)
    if "-nStep" in sys.argv:
        i = sys.argv.index("-nStep")
        nStep = int(sys.argv[i+1])
        sys.argv.pop(i)
        sys.argv.pop(i)

    inputFile = sys.argv[1]
    criteria = sys.argv[3]
    configIterator = parserGens.parseEnsemble(inputFile+".rmsd")
    rmsdTime = array([rmsd for rmsd in configIterator][-1])

    if scale == None:
        if "OUTCAR" in inputFile:
            scale = 5E-15 #5fs (potim=5)
        else:
            scale = 0.001E-12 #0.001ps (default for units = metal)

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
        steps,isfs = ISFFullRmsd(atomsTime,basis,rmsdTime,criteria=criteria,nqVecs=3,nStep=nStep)
        ylab = "$F_T(q,t)$"
        pl.title("Intermediate Scattering Function RMSD"+criteria)


    elif isfType == "self":
        steps,isfs = ISFSelfRmsd(atomsTime,basis,rmsdTime,criteria=criteria,nqVecs=3,nStep=nStep)
        if logtEnable:
            steps=log10(steps)
        pl.plot(steps,isfs)
        ylab = "$F_S(q,t)$"

    scale /= 1E-12 #picosecond conversion
    steps = [i*scale for i in steps]
    if logtEnable:
        steps=log10(steps)
    
    if logtEnable:
        header = "#Generated by %s using file %s\n log(t) ISF \n"%(sys.argv[0].split("/")[-1],inputFile)
        pl.xlabel("$log_{10}(t (ps))$")
    else:
        header = "#Generated by %s using file %s\n t(ps) ISF \n"%(sys.argv[0].split("/")[-1],inputFile)
        pl.xlabel("t(ps)")
    pl.ylabel(ylab)
    
    criteria = criteria.replace("gt",">")
    criteria = criteria.replace("lt","<")
    outputFile = inputFile + ".isf" + isfType[0].upper() +"_"+criteria
    odata = header
    for b in range(len(steps)):
        odata += str(steps[b])+" "+str(isfs[b])+"\n"
    open(outputFile,"w").write(odata)

    if plotEnable:
        pl.show()
