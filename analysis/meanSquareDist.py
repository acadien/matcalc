#!/usr/bin/python

from scipy import weave,array,zeros
from scipy.weave import converters
from numpy import *

msdCode = """
double c,d;
double *seta,*setb;
for(int dt=1;dt<=Ntime/2;dt++){
  c=0.0;
  for(int i=0;i<Ntime-dt;i++){
    seta=&(atoms[i*Natom*3]);
    setb=&(atoms[(i+dt)*Natom*3]);
    for(int j=0;j<Natom;j++){     //loop over each atom
      for(int dof=0;dof<3;dof++){       //loop over each DOF
        d = seta[j*3+dof]-setb[j*3+dof];
        if(d>lengths[dof]/2.0) d -= lengths[dof];
        if(d<-lengths[dof]/2.0) d += lengths[dof];
        c = c + d*d;
      }
    }
  }
  msd[dt-1]=c/Natom/(Ntime-dt);
}
"""
#Calculates the pair wise dist between for 2 sets points enforcing periodic bounds, puts the results in the ds array
def meanSquareDist(atoms,Natom,Ntime,lengths):
    print atoms[0][5]
    atoms=atoms.ravel() #atoms is an Ntime x Natom array... now flattened
    print atoms[15],atoms
    delT=array(range(1,Ntime/2+1))
    msd=zeros(len(delT))
    weave.inline(msdCode,['atoms','Ntime','Natom','lengths','msd','delT'])
    return delT,msd
