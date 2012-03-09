#!/usr/bin/python
import sys
import pylab as pl
#mine
from diffractools import *

#What does this script do?
#Generates the diffraction pattern for a hexagonal system

def usage():
    print "%s <lambda> <a=b> <c> <min2theta> <max2theta>"%sys.argv[0]
    print "Assumes a perfect hexagonal lattice with alpha=beta=90, gamma=120"
    

if len(sys.argv)!=6:
    usage()
    exit()

#Lattice and XRay parameters
lamda=float(sys.argv[1])

a=b=float(sys.argv[2])
c=float(sys.argv[3])

alpha=beta=90
gamma=120

min2theta=float(sys.argv[4])
max2theta=float(sys.argv[5])

#Reciprical Space calculations
A1=a
A2=a+2*b
A3=c
dspaces=[calc_dspace(A1,A2,A3,alpha,beta,gamma,hkl) for hkl in hkls]
thetas=[2.0*calc_theta(lamda,d) for d in dspaces]

#thetas=[i for i in thetas if i<max2theta and i>min2theta]

[pl.plot([i,i],[0.0,5.0]) for i in thetas]
pl.xlim([min2theta,max2theta])
pl.show()
