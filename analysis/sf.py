

#!/usr/bin/python

from scipy import fftpack
from math import *
import cmath

#Calculate the Structure factor from the rdf
def sf(rdfX,rdfY,Lmax=20.0,qbins=1024,ndens=0.05):
    minq,maxq,dq=0.2,20.0,(20.0-0.2)/qbins
    qs=[i*dq+minq for i in range(qbins)]
    print rdfX
    #sf=[1+ndens*sum([rdfY[ir]*cmath.exp(complex(0.,q*r)) for ir,r in enumerate(rdfX)]) for i,q in enumerate(qs)]#/2pi?
    sf=[1+4*pi*ndens*sum([(rdfY[ir]-1.0)*sin(q*r)/(q*r) for ir,r in enumerate(rdfX)]) for i,q in enumerate(qs)]#/2pi?
        
    return qs,sf

def sf2(atoms,Lmax=20.0, qbins=1024, ndens=0.05):
    
