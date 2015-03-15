#!/usr/bin/python

import sys
import numpy as np
import scipy.optimize as optimize
import pylab as pl
#mine
import utils

#Cite PRE 52,4 october 1995, W. Kob and H. Anderson 

utils.usage(["<ISFs file> <cutoff>"],1,3)

def parseISF(isfFile):
    dat = open(isfFile,"r").readlines()[1:]
    head = dat.pop(0)
    time,isfs = zip(*map(lambda x: map(float,x.split()),dat))
    return time,isfs

def kww(t,beta,a,tao):
    return a*np.exp(- (t/tao)**beta)

time,isfs = parseISF(sys.argv[1])
cutoff = 1.0
if len(sys.argv)==3:
    cutoff=float(sys.argv[2])


p2 = sys.argv[3]
a1,b1 = parseISF(p2)
icut = sum([1 for i in time if i<cutoff])
print a1,b1

#Fitting
beta,a,tao = optimize.curve_fit(kww,time[icut:],isfs[icut:])[0]
print "Beta = %f\n A = %f\n tao = %f\n"%(beta,a,tao)

#Plotting
#pl.plot(time[:icut],isfs[:icut])
pl.scatter(time,isfs,marker="x",lw=2,c="black",alpha=0.8)
longTime = range(2,10**5)
pl.plot(longTime,kww(longTime,beta,a,tao))

pl.plot(a1,b1,c="red")

pl.show()
