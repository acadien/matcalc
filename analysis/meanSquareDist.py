#!/usr/bin/python

from scipy import weave
from scipy.weave import converters
from numpy import *

msdCode = """
double c,d,atomd;
double *seta,*setb;
int dt;

#pragma omp parallel for default(shared) private(dt,c)
for(dt=1;dt<Ntime;dt++){

  c=0.0;

  for(int i=0;i<Ntime-dt;i++){

    seta=&(atoms[i*Natom*3]);
    setb=&(atoms[(i+dt)*Natom*3]);

    for(int j=0;j<Natom;j++){
      atomd=0
      for(int dof=0;dof<3;dof++){       //loop over each Degree of Freedom

        d = seta[j*3+dof]-setb[j*3+dof];
        if( d > lengths[dof]/2.0 ) d -= lengths[dof];  //Lower PBC
        if( d <-lengths[dof]/2.0 ) d += lengths[dof];  //Upper PBC

        c += d*d;
        atomd += d*d;
      }
      msdByAtom[(dt-1)*Natom + j] = atomd
    }
  }

  msd[dt-1] = c/(Natom*(Ntime-dt));

}
"""

#Calculates the MSD, doesn't take into big atomic jumps due to PBC
def meanSquareDist(atoms,Natom,Ntime,lengths,byAtom=False):
    atoms=array(atoms).ravel() #atoms is an Ntime x Natom x 3 array... now flattened
    delT=array(range(1,Ntime+1))
    msd=zeros(len(delT))
    msdByAtom=zeros(len(delT)*Natom)

    compiler_args=['-march=native -O3 -fopenmp']
    headers=r"""#include <omp.h>"""
    libs=['gomp']

    weave.inline(msdCode,['atoms','Ntime','Natom','lengths','msd','msdByAtom'],\
                     extra_compile_args=compiler_args,\
                     support_code=headers,\
                     libraries=libs)

    msdByAtom.shape=[len(delT),Natom]

    if byAtom:
        return delt,msd,msdByAtom
    return delT,msd

undoPBCcode = """
int index;
double d;
double *seta,*setb;

#pragma omp parallel for default(shared) private(index,d)
//first you need to unwrap the periodic boundary conditions in certain cases to ensure continuity
for(int i=0;i<Ntime-1;i++){

  seta=&(atoms[i*Natom*3]);
  setb=&(atoms[(i+1)*Natom*3]);

  for(int j=0;j<Natom;j++){
    for(int dof=0;dof<3;dof++){       //loop over each Degree of Freedom
      index=j*3+dof;
      d = seta[index]-setb[index];
      if( d > lengths[dof]/2.0 ) setb[index] += lengths[dof];  //Lower PBC
      if( d <-lengths[dof]/2.0 ) setb[index] -= lengths[dof];  //Upper PBC
    }
  }
} 

"""

msdRefCode = """
double c,d,atomd;
double *seta,*setb;

seta=&(refAtoms[0]);

//#pragma omp parallel for default(shared) private(c,d,atomd)
for(int i=0;i<Ntime;i++){

  setb=&(atoms[(i)*Natom*3]);

  c=0.0;
  for(int j=0;j<Natom;j++){
    atomd = 0;
    for(int dof=0;dof<3;dof++){       //loop over each Degree of Freedom

      d = seta[j*3+dof]-setb[j*3+dof];
      //if( d > lengths[dof]/2.0 ) d -= lengths[dof];  //Lower PBC
      //if( d <-lengths[dof]/2.0 ) d += lengths[dof];  //Upper PBC

      atomd += d*d;
    }
    c += atomd;
    msdByAtom[i*Natom + j] = atomd;
  }

  msd[i] = c/(Natom);//*Ntime);

}
"""

#Calculates the MSD enforcing periodic boundary conditions
#ref is an index into the atoms array
def meanSquareDistRef(atoms,ref,Natom,Ntime,lengths,byAtom=False):

    atoms=array(atoms).ravel()       #atoms is an Ntime x Natom x 3 array... now flattened
    delT=array(range(1,Ntime+1))
    msd=zeros(len(delT))
    msdByAtom=zeros(Ntime*Natom)

    compiler_args=['-march=native -O3 -fopenmp']
    headers=r"""#include <omp.h>"""
    libs=['gomp']

    weave.inline(undoPBCcode,['atoms','Ntime','Natom','lengths'],\
                     extra_compile_args=compiler_args,\
                     support_code=headers,\
                     libraries=libs)

    #do this after undoing PBC to ensure refAtoms are in same state as atoms set.
    atoms.shape=(Ntime,Natom,3)
    refAtoms=atoms[ref] #refAtoms is an Natom x 3 array... now flattened
    atom=atoms.ravel()
    refAtoms=refAtoms.ravel()
    weave.inline(msdRefCode,['atoms','refAtoms','Ntime','Natom','lengths','msd','msdByAtom'],\
                     extra_compile_args=compiler_args,\
                     support_code=headers,\
                     libraries=libs)
    msdByAtom.shape=[Ntime,Natom]

    if byAtom:
        return delT,msd,msdByAtom
    return delT,msd


msdRefAtomCode = """
double c,d;
double *seta,*setb;

seta=&(refAtoms[0]);

#pragma omp parallel for default(shared) private(c,d)
for(int i=0;i<Ntime;i++){

  setb=&(atoms[(i)*Natom*3]);

  for(int j=0;j<Natom;j++){
    c=0.0;
    for(int dof=0;dof<3;dof++){       //loop over each Degree of Freedom

      d = seta[j*3+dof]-setb[j*3+dof];
      //if( d > lengths[dof]/2.0 ) d -= lengths[dof];  //Lower PBC
      //if( d <-lengths[dof]/2.0 ) d += lengths[dof];  //Upper PBC

      c += d*d;

    }

    msdAtom[i*Natom+j] = c;

  }

}
"""

#Calculates the MSD enforcing periodic boundary conditions, leaves info in terms of per atom
def meanSquareDistRefAtom(atoms,ref,Natom,Ntime,lengths):

    atoms=array(atoms).ravel() #atoms is an Ntime x Natom x 3 array... now flattened
    delT=array(range(1,Ntime+1))
    msdAtom=zeros(Ntime*Natom)

    compiler_args=['-march=native -O3 -fopenmp']
    headers=r"""#include <omp.h>"""
    libs=['gomp']

    weave.inline(undoPBCcode,['atoms','Ntime','Natom','lengths'],\
                     extra_compile_args=compiler_args,\
                     support_code=headers,\
                     libraries=libs)

    #do this after undoing PBC to ensure refAtoms are in same state as atoms set.
    atoms.shape=(Ntime,Natom,3)
    refAtoms=atoms[ref] #refAtoms is an Natom x 3 array... now flattened
    atom=atoms.ravel()
    refAtoms=refAtoms.ravel()
    weave.inline(msdRefAtomCode,['atoms','refAtoms','Ntime','Natom','lengths','msdAtom'],\
                     extra_compile_args=compiler_args,\
                     support_code=headers,\
                     libraries=libs)
    msdAtom.shape=[Ntime,Natom]
    msdAtom=msdAtom.T
    return delT,msdAtom
