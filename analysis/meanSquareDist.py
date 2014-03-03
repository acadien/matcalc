#!/usr/bin/python

from scipy import weave
from scipy.weave import converters
from numpy import *

#Yeah that's right stick your code in a fucking comment cause its a weave hack.
msdCode = """
double c,d;
double *seta,*setb;
int dt;

#pragma omp parallel for default(shared) private(dt,c)
for(dt=1;dt<Ntime;dt++){

  c=0.0;

  for(int i=0;i<Ntime-dt;i++){

    seta=&(atoms[i*Natom*3]);
    setb=&(atoms[(i+dt)*Natom*3]);

    for(int j=0;j<Natom;j++){
      for(int dof=0;dof<3;dof++){       //loop over each Degree of Freedom

        d = seta[j*3+dof]-setb[j*3+dof];
        if( d > lengths[dof]/2.0 ) d -= lengths[dof];  //Lower PBC
        if( d <-lengths[dof]/2.0 ) d += lengths[dof];  //Upper PBC

        c += d*d;

      }
    }
  }

  msd[dt-1] = c/(Natom*(Ntime-dt));

}
"""

#Calculates the MSD enforcing periodic boundary conditions
def meanSquareDist(atoms,Natom,Ntime,lengths):
    atoms=array(atoms).ravel() #atoms is an Ntime x Natom x 3 array... now flattened
    delT=array(range(1,Ntime+1))
    msd=zeros(len(delT))

    compiler_args=['-march=native -O3 -fopenmp']
    headers=r"""#include <omp.h>"""
    libs=['gomp']

    weave.inline(msdCode,['atoms','Ntime','Natom','lengths','msd'],\
                     extra_compile_args=compiler_args,\
                     support_code=headers,\
                     libraries=libs)
    return delT,msd
