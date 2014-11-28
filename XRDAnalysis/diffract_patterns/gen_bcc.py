#!/usr/bin/python
import sys
import pylab as pl
#mine
from diffractools import *

#What does this script do?
#Generates the diffraction pattern for a hexagonal system

def usage():
    print "%s <lambda> <a> <max2theta>"%sys.argv[0]
    print "Assumes a perfect bcc lattice with a=b=c and alpha=beta=gamma=90"
    

if len(sys.argv)!=4:
    usage()
    exit()

#Lattice and XRay parameters
lamda=float(sys.argv[1])

a=b=c=float(sys.argv[2])

alpha=beta=gamma=90

max2theta=float(sys.argv[3])

#Reciprical Space calculations
dspaces=[calc_dspace(a,b,c,alpha,beta,gamma,hkl) for hkl in hkls]
thetas=[calc_theta(lamda,d) for d in dspaces]
thetas=[i for i in thetas if i<max2theta and i>0]

[pl.plot([2*i,2*i],[0.0,5.0]) for i in thetas]
pl.show()
