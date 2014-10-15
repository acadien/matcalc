#!/usr/bin/python

from scipy import weave
from scipy.weave import converters
from numpy import *


#############################################################################################
#############################################################################################
#############################################################################################

undoPBCcode = """
int index;
double d;
double *seta,*setb;

//first you need to unwrap the periodic boundary conditions in certain cases to ensure continuity
for(int i=0;i<Ntime-1;i++){

  seta=&(atoms[i*Natom*3]);
  setb=&(atoms[(i+1)*Natom*3]);

  for(int j=0;j<Natom;j++){
    for(int dof=0;dof<3;dof++){       //loop over each Degree of Freedom
      index=j*3+dof;
      d = seta[index]-setb[index];
      if( d > lengths[dof]/3.0 ) setb[index] += lengths[dof];  //Lower PBC
      if( d <-lengths[dof]/3.0 ) setb[index] -= lengths[dof];  //Upper PBC
    }
  }
} 

"""

rmsdRefCode = """
double c,d,atomd;
double *seta,*setb;

seta=&(refAtoms[0]);

for(int i=0;i<Ntime;i++){

  setb=&(atoms[(i)*(int)Natom*3]);

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
    rmsdByAtom[i*Natom + j] = sqrt(atomd);
  }

  rmsd[i] = sqrt(c/(Natom));

}
"""

#Calculates the RMSD enforcing periodic boundary conditions
#ref is an index into the atoms array
def rootMeanSquareDistRef(atoms,ref,Natom,Ntime,lengths,byAtom=False):
    Natom=int(Natom)
    atoms=array(atoms).ravel()       #atoms is an Ntime x Natom x 3 array... now flattened
    delT=array(range(1,Ntime+1))
    rmsd=zeros(len(delT))
    rmsdByAtom=zeros(Ntime*Natom)

#    compiler_args=['-march=native -O3 -fopenmp']
#    headers=r"""#include <omp.h>"""
#    libs=['math.h']

    weave.inline(undoPBCcode,['atoms','Ntime','Natom','lengths'])#,libraries=libs)
#                     extra_compile_args=compiler_args,\
#                     support_code=headers,\


    #do this after undoing PBC to ensure refAtoms are in same state as atoms set.
    atoms.shape=(Ntime,Natom,3)
    refAtoms=atoms[ref] #refAtoms is an Natom x 3 array... now flattened
    atom=atoms.ravel()
    refAtoms=refAtoms.ravel()
    weave.inline(rmsdRefCode,['atoms','refAtoms','Ntime','Natom','lengths','rmsd','rmsdByAtom'])#,libraries=libs)
#                     extra_compile_args=compiler_args,\
#                     support_code=headers,\

    rmsdByAtom.shape=[Ntime,Natom]

    if byAtom:
        return delT,rmsd,rmsdByAtom
    return delT,rmsd

#############################################################################################
#############################################################################################
#############################################################################################

rmsdWindowCode = """
double c,d,atomd;
double *seta,*setb;

for(int i=0;i<Ntime;i++){

  if(i<window)
    seta=&(atoms[0]);
  else
    seta=&(atoms[(i-window)*(int)Natom*3]);

  setb=&(atoms[(i)*(int)Natom*3]);

  c=0.0;
  for(int j=0;j<Natom;j++){
    atomd = 0;
    for(int dof=0;dof<3;dof++){       //loop over each Degree of Freedom

      d = seta[j*3+dof]-setb[j*3+dof];
      atomd += d*d;
    }
    c += atomd;
    rmsdByAtom[i*Natom + j] = sqrt(atomd);
  }

  rmsd[i] = sqrt(c/(Natom));

}

"""
#Calculates the RMSD enforcing periodic boundary conditions over a window size
#'window' is an integer that defines the seperation between seta and setb atoms when calculating RMSD.
def rootMeanSquareDistWindow(atoms,window,Natom,Ntime,lengths,byAtom=False):
    Natom=int(Natom)
    atoms=array(atoms).ravel()       #atoms is an Ntime x Natom x 3 array... now flattened
    delT=array(range(1,Ntime+1))
    rmsd=zeros(len(delT))
    rmsdByAtom=zeros(Ntime*Natom)

    weave.inline(undoPBCcode,['atoms','Ntime','Natom','lengths'])#,libraries=libs)

    #do this after undoing PBC to ensure refAtoms are in same state as atoms set.
    atoms.shape=(Ntime,Natom,3)
    atom=atoms.ravel()
    weave.inline(rmsdWindowCode,['atoms','window','Ntime','Natom','lengths','rmsd','rmsdByAtom'])#,libraries=libs)

    rmsdByAtom.shape=[Ntime,Natom]

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

    rmsdAtom[i*Natom+j] = sqrt(c);

  }

}
"""

#Calculates the RMSD enforcing periodic boundary conditions, leaves info in terms of per atom
def rootMeanSquareDistRefAtom(atoms,ref,Natom,Ntime,lengths):

    atoms=array(atoms).ravel() #atoms is an Ntime x Natom x 3 array... now flattened
    delT=array(range(1,Ntime+1))
    rmsdAtom=zeros(Ntime*Natom)

#    compiler_args=['-march=native -O3 -fopenmp']
#    headers=r"""#include <omp.h>"""
#    libs=['math.h']

    weave.inline(undoPBCcode,['atoms','Ntime','Natom','lengths'])#,libraries=libs)
#                     extra_compile_args=compiler_args,\
#                     support_code=headers,\

    #do this after undoing PBC to ensure refAtoms are in same state as atoms set.
    atoms.shape=(Ntime,Natom,3)
    refAtoms=atoms[ref] #refAtoms is an Natom x 3 array... now flattened
    atom=atoms.ravel()
    refAtoms=refAtoms.ravel()
    weave.inline(rmsdRefAtomCode,['atoms','refAtoms','Ntime','Natom','lengths','rmsdAtom'])#,libraries=libs)
#                     extra_compile_args=compiler_args,\
#                     support_code=headers,\

    rmsdAtom.shape=[Ntime,Natom]
    rmsdAtom=rmsdAtom.T
    return delT,rmsdAtom

#############################################################################################
#############################################################################################
#############################################################################################

def unwrap(atoms,Natom,Ntime,lengths):
    atoms=array(atoms).ravel() #atoms is an Ntime x Natom x 3 array... now flattened
    rmsdAtom=zeros(Ntime*Natom)
    weave.inline(undoPBCcode,['atoms','Ntime','Natom','lengths'])
    atoms.shape=(Ntime,Natom,3)
    return atoms




#garbage below



"""
rmsdCode = ""
double c,d,atomd;
double *seta,*setb;
int dt;

for(dt=1;dt<Ntime;dt++){

  c=0.0;

  for(int i=0;i<Ntime-dt;i++){

    seta=&(atoms[i*Natom*3]);
    setb=&(atoms[(i+dt)*Natom*3]);

    for(int j=0;j<(int)Natom;j++){
      atomd=0
      for(int dof=0;dof<3;dof++){       //loop over each Degree of Freedom

        d = seta[j*3+dof]-setb[j*3+dof];
        if( d > lengths[dof]/2.0 ) d -= lengths[dof];  //Lower PBC
        if( d <-lengths[dof]/2.0 ) d += lengths[dof];  //Upper PBC

        c += d*d;
        atomd += d*d;
      }
      rmsdByAtom[(dt-1)*Natom + j] = sqrt(atomd)
    }
  }

  rmsd[dt-1] = sqrt(c/(Natom*(Ntime-dt)));

}
""

#Calculates the RMSD, doesn't take into big atomic jumps due to PBC
def rootMeanSquareDist(atoms,Natom,Ntime,lengths,byAtom=False):
    Natom=int(Natom)
    atoms=array(atoms).ravel() #atoms is an Ntime x Natom x 3 array... now flattened
    delT=array(range(1,Ntime+1))
    rmsd=zeros(len(delT))
    rmsdByAtom=zeros(len(delT)*Natom)

#    compiler_args=['-march=native -O3 -fopenmp']
#    headers=r""#include <omp.h>""
#    libs=['gomp']

    weave.inline(rmsdCode,['atoms','Ntime','Natom','lengths','rmsd','rmsdByAtom'])
#                     extra_compile_args=compiler_args,\
#                     support_code=headers,\
#                     libraries=libs)

    rmsdByAtom.shape=[len(delT),Natom]

    if byAtom:
        return delt,rmsd,rmsdByAtom
    return delT,rmsd
"""
