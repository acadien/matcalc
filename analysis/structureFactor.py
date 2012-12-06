#!/usr/bin/python

from scipy import fftpack
from math import *
import cmath
#mine
#import paircor


def structureFactor(pairCorX,pairCorY,Lmax=20.0,qbins=1024,ndens=0.05):
    minq,maxq,dq=0.2,20.0,(20.0-0.2)/qbins
    qs=[i*dq+minq for i in range(qbins)]
    print pairCorX
    sf=[1+ndens*sum([pairCorY[ir]*cmath.exp(complex(0.,q*r)) for ir,r in enumerate(pairCorX)]) for i,q in enumerate(qs)]#/2pi?
        

    return qs,sf
    
