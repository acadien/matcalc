#!/usr/bin/python

import sys
import numpy as np
import scipy.optimize as optimize
import pylab as pl
#mine
import utils

utils.usage(["<ISFs file> <cutoff>"],1,3)

def parseISF(isfFile):
    dat = open(isfFile,"r").readlines()[1:]
    head = dat.pop(0)
    time,isfs = zip(*map(lambda x: map(float,x.split()),dat))
    return time,isfs

def vtf(t,beta,a,tao):
    return a*np.exp(- (t/tao)**beta)

time,isfs = parseISF(sys.argv[1])
cutoff = 1.0
if len(sys.argv)==3:
    cutoff=float(sys.argv[2])

icut = sum([1 for i in time if i<cutoff])

#Fitting
beta,a,tao = optimize.curve_fit(kww,time[icut:],isfs[icut:])[0]
print "Beta = %f\n A = %f\n tao = %f\n"%(beta,a,tao)

#Plotting
#pl.plot(time[:icut],isfs[:icut])
pl.scatter(time,isfs,marker="x",lw=2,c="black",alpha=0.8)
#longTime = [i/100. for i in range(cut*100,10**6)]
longTime = list(time)+range(int(time[-1]),int(time[-1])*100)
pl.plot(longTime,kww(longTime,beta,a,tao))

data=map(lambda x: " ".join(map(str,x))+"\n",zip(longTime,kww(longTime,beta,a,tao)))
open(sys.argv[1]+".kww","w").writelines(data)

pl.show()
