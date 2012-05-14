#!/usr/bin/python
import sys
import pylab as pl
import re
import os
import numpy as np
from numpy import *
from math import asin

def readchi(chifile):
    i=0
    for line in chifile:
        if i<3:
            i+=1
            continue
        line=line.lstrip()
        line=line.rstrip()
        if i==3:
            i+=1
            length=int(line)
            break
    xxs=np.zeros(length)
    yys=np.zeros(length)
    i=0;
    for line in chifile:
        line=line.strip().split()
        xxs[i]=line[0]
        yys[i]=line[1]
        i+=1
    return xxs,yys

def theta2q(theta2,lamda):
    return [4.0*pi*sin(radians(i/2.0))/lamda for i in theta2] 

#Return 2theta in degrees
def q2theta(q,lamda):
    return degrees([2*asin(i*lamda/4.0/pi) for i in q])


if __name__=="__main__":
    if len(sys.argv) < 3:
        print("Bad arguments, proper usage:")
        print("./chi2q.py <chifile1..N> <lambda>")
        print "Converts chi file to q file which is intensity vs q value"
        exit(0)

    lamda=float(sys.argv[-1])
    chifilenames=sys.argv[1:-1]
    qfilenames=[i.split(".")[0]+".q" for i in chifilenames]

    #####################
    #Read in the chi file
    for i,qfilename in enumerate(qfilenames):
        x,y=readchi(open(chifilenames[i],"r"))
        x=theta2q(x,lamda)

        qfil=open(qfilename,"w")
        qfil.write("Q in Angstom^-1 vs Intensity\n")
        qfil.write("%d\n"%len(x))
        for a,b in zip(x,y):
            qfil.write("%5.5f %5.5f\n"%(a,b))
        qfil.close()


        
    
