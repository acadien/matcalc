#!/usr/bin/python

import sys
import pylab as pl
from scipy import *

def usage():
    print "Usage: %s <PCDAT file> <time divisions (10)> <start iter(0)> <end iter(-1)>"%sys.argv[0]

starti=0
endi=None
divs=10
enableAvg=False
if len(sys.argv)<2:
    usage()
    exit(0)
if len(sys.argv)>=3:
    divs=int(sys.argv[2])
if len(sys.argv)>=4:
    starti=int(sys.argv[3])
if len(sys.argv)>=5:
    if sys.argv[4]=='-1':
        endi=None
    else:
        endi=int(sys.argv[4])

pcfile=sys.argv[1]
pcdat=open(pcfile,"r").readlines()
head=pcdat[:12]
Nbins=int(pcdat[6])
delr=float(pcdat[8])*1E10

pcdat=pcdat[11:]

xx=[i*delr for i in range(Nbins)]

yy=list()

if endi!=None:
    pcdat=pcdat[(Nbins+1)*starti:(Nbins+1)*(endi+1)]
else:
    pcdat=pcdat[(Nbins+1)*starti:]

#Read in bin info
yy=[map(float,pcdat[i*Nbins+i+1:(i+1)*Nbins+i+1]) for i in range(len(pcdat)/(Nbins+1))]
Ntot=len(yy)
Nspace=Ntot/divs

offset=max([max(i) for i in yy])/20.

if len(yy)!=1:
    starts=[i*Nspace for i in range(divs)]
    ends=[min([i+Nspace,Ntot]) for i in starts]
    labs=["[%d,%d]"%(a+starti,b+starti) for a,b in zip(starts,ends)]
    for i in range(divs):
        si=starts[i]
        ei=ends[i]
        ydiv = [sum(y)/len(y)+offset*i for y in zip(*yy[si:ei])]
        pl.plot(xx,ydiv,label=labs[i])
else:
    labs="1"
    pl.plot(xx,yy[0],label=labs)
pl.xlabel("r ($\AA$)")
pl.ylabel("g(r)")
if endi==None:
    endi=Ntot
pl.title("%s [%d,%d]"%(pcfile,starti,endi))
pl.legend(labs,title="Iterations",loc=0)
pl.show()

