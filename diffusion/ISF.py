#!/usr/bin/python

#Takes a dump.dat file, parses it, generates van-hove, fourier transform to
#generate intermediate scattering function

#SHOULD BE RUN ONLY ON NVT/NVE SIMULATIONS WITH UNWRAPPED COORDINATES
#Unwraps outcar coordinates

#mine
import plotRemote as pr
#theirs
import sys
import itertools,random
from scipy import weave
from scipy.weave import converters
from scipy.integrate import simps
from numpy import array,zeros,pi,log10,cos,sin,sqrt
import pylab as pl
#mine
import utils
import parserGens,lammpsIO,outcarIO
from rootMeanSquareDist import unwrap
from vanHove import vanHoveSelf


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

    qxs[qlen] = (2*a*sqrt(1-a2-b2))*q;
    qys[qlen] = (2*b*sqrt(1-a2-b2))*q;
    qzs[qlen] = (1-2*(a2+b2))*q;

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

def ISFSelf(atoms,basis,steps=None,nqVecs=1,qmax=1.9):#qmax (first peak=1.9, second peak=3.32)
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

    nStep = steps.shape[0]

    isfs = zeros(nStep)
    q=qmax

    weave.inline(ISFSelfSphereRefCode,['nqVecs','steps','nStep','isfs','atoms','nTime','nAtom','q'])
    return steps,isfs


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
while( qlen++ < nqVecs ){
    a = rand()*2/(double)RAND_MAX-1.0;
    b = rand()*2/(double)RAND_MAX-1.0;
    a2 = a*a;
    b2 = b*b;

    if(a2 + b2 > 1.0) 
        continue;

    qx = (2*a*sqrt(1-a2-b2))*q;
    qy = (2*b*sqrt(1-a2-b2))*q;
    qz = (1-2*(a2+b2))*q;

    for(int s=0; s<nStep; s++){
      for(int t=0; t<nTime-stepSize; t++){
        seta = &(atoms[(int)(t*nAtom*3)]);
        setb = &(atoms[(int)((t+stepSize)*nAtom*3)]);

        for(int i=0; i<nAtom; i++){
          da = seta[i*3+0]*qx + seta[i*3+1]*qy + seta[i*3+2]*qz;
          db = setb[i*3+0]*qx + setb[i*3+1]*qy + setb[i*3+2]*qz;

          isfs[s] += cos(da)*cos(db);// + sin(da)*sin(db);
    }}}
}

/*
for(int t=0; t<nTime; t++){
  rhoqcos[t] /= nqVecs;
  rhoqsin[t] /= nqVecs;
}

for(int s=0; s<nStep; s++){ //loop over step size
  stepSize = steps[s];
  
  for(int t=0; t<nTime-stepSize; t++){ //loop over time steps
    isfs[s] += rhoqcos[t];//(rhoqcos[t+stepSize]*rhoqcos[t] + rhoqsin[t+stepSize]*rhoqsin[t]) / (nTime-stepSize);
  }
}*/
"""

def ISFFull(atoms,basis,steps=None,nqVecs=1,qmax=1.9):#kmax (first peak=1.9, second peak=3.32)
    atoms = array(atoms)
    basis = array(basis)
    nTime = atoms.shape[0]
    nAtom = atoms.shape[1]

    if steps == None:
        nSteps = 200
        steps = sorted(list(set(array([int(10**(t*10.0/nSteps-1.0)) for t in range(nSteps)]))))
        steps.remove(0)
        steps = [t for t in steps if t<nTime]
        steps=array(steps)

    atoms = atoms.ravel()

    nStep = steps.shape[0]

    isfs = zeros(nStep)
    q=qmax

    weave.inline(ISFFullSphereRefCode,['nqVecs','steps','nStep','isfs','atoms','nTime','nAtom','q'])
    return steps,isfs


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
        steps,isfs = ISFFull(atomsTime,basis,nqVecs=3)
        """
        q=1.9
        pq=list()
        for a1 in atomsTime:
            pqn=0
            for n in range(1):
                r=1.1
                while(r>1.0):
                    a = random.random()*2-1.0
                    b = random.random()*2-1.0
                    r = a*a+b*b
                qx = (2*a*sqrt(1-a*a-b*b))*sqrt(q)
                qy = (2*b*sqrt(1-a*a-b*b))*sqrt(q)
                qz = (1-2*(a*a+b*b))*sqrt(q)
                d = [(qx*a[0]+qy*a[1]+qz*a[2]) for a in a1]
                pqn += sum([sin(a)**2+cos(a)**2 for a in d])
            pq.append(pqn/1.0)
          pl.plot(pq)
          pl.show()
          """
        import pylab as pl
        print isfs
        pl.plot(steps,isfs)
        pl.show()
        exit(0)
            
        #rbins,bins = vanHoveTotal(atomsTime,basis,steps,cutr=cutr,nBin=nBin,norm=norm)
        #ylab = "$G(r,t) / V_{shell}(r)$"
        pass
    elif isfType == "self":
        steps,isfs = ISFSelf(atomsTime,basis,nqVecs=3)
        if logtEnable:
            steps=log10(steps)
        pl.plot(steps,isfs)
        if logtEnable:
            header = "#Generated by %s using file %s\n log(t) ISF \n"%(sys.argv[0].split("/")[-1],inputFile)
            pl.xlabel("$log_{10}(t)$")
        else:
            header = "#Generated by %s using file %s\n t ISF \n"%(sys.argv[0].split("/")[-1],inputFile)
            pl.xlabel("t")
        pl.ylabel("$F_S(q,t)$")
    
    outputFile = inputFile + ".isf" + isfType[0].upper()
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
