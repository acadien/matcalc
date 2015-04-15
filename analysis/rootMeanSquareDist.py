#!/usr/bin/python

from scipy import weave
from scipy.weave import converters
from numpy import *


#############################################################################################
#############################################################################################
#############################################################################################

undoPBCcode = """
int xIndex,yIndex,zIndex;
int minDof1,minDof2,minDof3;
double dx,dy,dz,c,minC;
double *seta,*setb;

//first you need to unwrap the periodic boundary conditions in certain cases to ensure continuity
for(int i=0;i<nTime-1;i++){
  seta=&(atoms[(i  )*nAtom*3]);
  setb=&(atoms[(i+1)*nAtom*3]);

  for(int j=0;j<nAtom;j++){

    xIndex=j*3+0;
    yIndex=j*3+1;
    zIndex=j*3+2;

    //Find the minimum image distance between atom and itself in next timestep.
    minC = 100000.0;
    for(int dof1=-1;dof1<2;dof1++){       
    for(int dof2=-1;dof2<2;dof2++){       
    for(int dof3=-1;dof3<2;dof3++){       
      dx = seta[xIndex] - (setb[xIndex] + dof1*b[0] + dof2*b[3] + dof3*b[6]);
      dy = seta[yIndex] - (setb[yIndex] + dof1*b[1] + dof2*b[4] + dof3*b[7]);
      dz = seta[zIndex] - (setb[zIndex] + dof1*b[2] + dof2*b[5] + dof3*b[8]);
      c = dx*dx+dy*dy+dz*dz;

      if(c<minC){
        minC = c;
        minDof1 = dof1;
        minDof2 = dof2;
        minDof3 = dof3;
      }
    }}}

    //If unwrapping is necessary apply it to all future steps for the atom.
    if(minDof1 != 0 || minDof2 != 0 || minDof3 != 0){
      dx = minDof1*b[0] + minDof2*b[3] + minDof3*b[6];
      dy = minDof1*b[1] + minDof2*b[4] + minDof3*b[7];
      dz = minDof1*b[2] + minDof2*b[5] + minDof3*b[8];

      for(int k=i+1; k<nTime;k++){
        atoms[(k*nAtom+j)*3]   += dx;
        atoms[(k*nAtom+j)*3+1] += dy;
        atoms[(k*nAtom+j)*3+2] += dz;
      }
    }
  }
}
"""

rmsdRefCode = """
double c,d,atomd;
double *seta,*setb;

seta=&(refAtoms[0]);

for(int i=0;i<nTime;i++){

  setb=&(atoms[(i)*(int)nAtom*3]);

  c=0.0;
  for(int j=0;j<nAtom;j++){
    atomd = 0;
    for(int dof=0;dof<3;dof++){       //loop over each Degree of Freedom

      d = seta[j*3+dof]-setb[j*3+dof];
      atomd += d*d;
    }
    c += atomd;
    rmsdByAtom[i*nAtom + j] = sqrt(atomd);
  }

  rmsd[i] = sqrt(c/(nAtom));

}
"""

#Calculates the RMSD enforcing periodic boundary conditions
#ref is an index into the atoms array
def rootMeanSquareDistRef(atoms,ref,basis,byAtom=False):
    atoms = array(atoms)
    basis = array(basis)

    nTime = atoms.shape[0]
    nAtom = atoms.shape[1]
    atoms=atoms.ravel()       #atoms is an nTime x nAtom x 3 array... now flattened
    b = basis.ravel()
    delT=linspace(1,nTime+1,nTime+1)
    rmsd=zeros(len(delT))
    rmsdByAtom=zeros(nTime*nAtom)
    debug = zeros(4)

#    compiler_args=['-march=native -O3 -fopenmp']
#    headers=r"""#include <omp.h>"""
#    libs=['math.h']
    weave.inline(undoPBCcode,['atoms','nTime','nAtom','b','debug'])#,libraries=libs)
#                     extra_compile_args=compiler_args,\
#                     support_code=headers,\

    #do this after undoing PBC to ensure refAtoms are in same state as atoms set.
    atoms.shape=(nTime,nAtom,3)
    refAtoms=atoms[ref] #refAtoms is an nAtom x 3 array... now flattened
    atom=atoms.ravel()
    refAtoms=refAtoms.ravel()
    weave.inline(rmsdRefCode,['atoms','refAtoms','nTime','nAtom','rmsd','rmsdByAtom'])#,libraries=libs)
#                     extra_compile_args=compiler_args,\
#                     support_code=headers,\

    rmsdByAtom.shape=[nTime,nAtom]

    if byAtom:
        return delT,rmsd,rmsdByAtom
    return delT,rmsd

#############################################################################################
#############################################################################################
#############################################################################################

rmsdWindowCode = """
double c,d,atomd;
double *seta,*setb;

for(int i=1;i<nTime;i++){

  if(i<window)
    seta=&(atoms[0]);
  else
    seta=&(atoms[(i-window)*(int)nAtom*3]);

  setb=&(atoms[(i)*(int)nAtom*3]);

  c=0.0;
  for(int j=0;j<nAtom;j++){
    atomd = 0;
    for(int dof=0;dof<3;dof++){       //loop over each Degree of Freedom
      d = seta[j*3+dof]-setb[j*3+dof];
      atomd += d*d;
    }
    c += atomd;
    rmsdByAtom[i*nAtom + j] = sqrt(atomd);
  }
  rmsd[i] = sqrt(c/(nAtom));

}

"""
#Calculates the RMSD enforcing periodic boundary conditions over a window size
#'window' is an integer that defines the seperation between seta and setb atoms when calculating RMSD.
def rootMeanSquareDistWindow(atoms,window,basis,byAtom=False):
    atoms = array(atoms)
    basis = array(basis)

    nTime = atoms.shape[0]
    nAtom = atoms.shape[1]
    atoms=atoms.ravel()       #atoms is an nTime x nAtom x 3 array... now flattened
    b = basis.ravel()
    delT=linspace(1,nTime+1,nTime+1)
    rmsd=zeros(len(delT))
    rmsdByAtom=zeros(nTime*nAtom)

    weave.inline(undoPBCcode,['atoms','nTime','nAtom','b'])#,libraries=libs)

    #do this after undoing PBC to ensure refAtoms are in same state as atoms set.
    atoms.shape=(nTime,nAtom,3)
    atom=atoms.ravel()
    weave.inline(rmsdWindowCode,['atoms','window','nTime','nAtom','rmsd','rmsdByAtom'])#,libraries=libs)

    rmsdByAtom.shape=[nTime,nAtom]

    if byAtom:
        return delT,rmsd,rmsdByAtom
    return delT,rmsd

#############################################################################################
#############################################################################################
#############################################################################################

rmsdRefAtomCode = """
double c,d;
double *seta,*setb;

seta=&(refAtoms[0]);

for(int i=0;i<nTime;i++){
  setb=&(atoms[(i)*nAtom*3]);

  for(int j=0;j<nAtom;j++){
    c=0.0;
    for(int dof=0;dof<3;dof++){       //loop over each Degree of Freedom
      d = seta[j*3+dof]-setb[j*3+dof];
      c += d*d;
    }
    rmsdAtom[i*nAtom+j] = sqrt(c);
  }
}
"""

#Calculates the RMSD enforcing periodic boundary conditions, leaves info in terms of per atom
def rootMeanSquareDistRefAtom(atoms,ref,basis):
    atoms = array(atoms)
    basis = array(basis)

    nTime = atoms.shape[0]
    nAtom = atoms.shape[1]
    atoms=atoms.ravel()       #atoms is an nTime x nAtom x 3 array... now flattened
    b = basis.ravel()
    delT=linspace(1,nTime+1,nTime+1)
    rmsd=zeros(len(delT))
    rmsdAtom=zeros(nTime*nAtom)
    debug = zeros(4)
#    compiler_args=['-march=native -O3 -fopenmp']
#    headers=r"""#include <omp.h>"""
#    libs=['math.h']

    weave.inline(undoPBCcode,['atoms','nTime','nAtom','b','debug'])#,libraries=libs)
#                     extra_compile_args=compiler_args,\
#                     support_code=headers,\
    #do this after undoing PBC to ensure refAtoms are in same state as atoms set.
    atoms.shape=(nTime,nAtom,3)
    refAtoms=atoms[ref] #refAtoms is an nAtom x 3 array... now flattened
    atom=atoms.ravel()
    refAtoms=refAtoms.ravel()
    weave.inline(rmsdRefAtomCode,['atoms','refAtoms','nTime','nAtom','rmsdAtom'])#,libraries=libs)
#                     extra_compile_args=compiler_args,\
#                     support_code=headers,\

    rmsdAtom.shape=[nTime,nAtom]
    rmsdAtom=rmsdAtom.T
    return delT,rmsdAtom

#############################################################################################
#############################################################################################
#############################################################################################

def unwrap(atoms,basis):

    atoms = array(atoms)
    nTime = atoms.shape[0]
    nAtom = atoms.shape[1]
    atoms = atoms.ravel() #atoms is an nTime x nAtom x 3 array... now flattened
    debug = zeros(4)

    b=array(basis).ravel()
    if(len(b)!=9): 
        print "noooooooooooooooo"
        exit(0)

    weave.inline(undoPBCcode,['atoms','nTime','nAtom','b','debug'])
    atoms.shape=(nTime,nAtom,3)
    return atoms

############################################################################################
#['atoms','nTime','nAtom','steps','nStep','rmsd']
rmsdCorrCode = """
double c,d,atomd;
double *seta,*setb;

int stepSize,k;

for(int s=0; s<nStep; s++){ //Loop over steps (rmsd[step] = )
  stepSize = steps[s];

  c = 0.0;
  k = 0;
  for(int i=0; i < nTime-stepSize; i+=stepSize){ //Loop over time by step size
  
    seta=&(atoms[(i         )*(int)nAtom*3]);
    setb=&(atoms[(i+stepSize)*(int)nAtom*3]);

    for(int j=0;j<nAtom;j++){
      atomd = 0;
      for(int dof=0;dof<3;dof++){       //loop over each Degree of Freedom
        d = seta[j*3+dof]-setb[j*3+dof];
        atomd += d*d;
      }
      c += sqrt(atomd);
      k++;
    }
  }
  rmsd[s] = c/k;
}
"""

#Calculates the RMSD enforcing periodic boundary conditions
#ref is an index into the atoms array
def rootMeanSquareCorrelate(atoms,basis,steps=None):
    atoms = array(atoms)
    basis = array(basis)
    nTime = atoms.shape[0]
    nAtom = atoms.shape[1]

    if steps==None:
        nStep=1000
        delTime = max( int(nTime/(nStep*10)) , 1)
        steps = array(range(1,nTime/10,delTime))
        nStep = steps.shape[0]
        print nStep
    else:
        nStep = len(steps)
        steps = array(steps)

    atoms=atoms.ravel()       #atoms is an nTime x nAtom x 3 array... now flattened
    b = basis.ravel()
    rmsd=zeros(nStep)
    debug = zeros(4)

    weave.inline(undoPBCcode,['atoms','nTime','nAtom','b','debug'])

    #do this after undoing PBC to ensure refAtoms are in same state as atoms set.
    weave.inline(rmsdCorrCode,['atoms','nTime','nAtom','steps','nStep','rmsd'])

    return steps,rmsd
