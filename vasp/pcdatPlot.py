#!/usr/bin/python

import sys
import pylab as pl
from scipy import *

def usage():
    print "Usage: %s <PCDAT file> <optional:start iter(0)> <optional:end iter(-1)>"%sys.argv[0]

starti=0
endi=None
if len(sys.argv)<1:
    usage()
    exit(0)
elif len(sys.argv)==3:
    starti=int(sys.argv[2])
elif len(sys.argv)==3:
    endi=int(sys.argv[3])

pcfile=sys.argv[1]
pcdat=open(pcfile,"r").readlines()
head=pcdat[:12]
Nbins=int(pcdat[6])
delr=float(pcdat[8])*1E10

pcdat=pcdat[11:]

xx=[i*delr for i in range(Nbins)]

yy=list()
while len(pcdat)>Nbins:
    pcdat.pop(0)
    yy.append(map(float,pcdat[0:Nbins]))
    pcdat=pcdat[Nbins:]

Ntot=len(yy)
Nused=len(yy[starti:endi])
yy = [sum(y)/len(y) for y in zip(*yy[starti:endi])]

print "Using %d iterations of %d possible"%(Nused,Ntot)


pl.plot(xx,yy)
pl.xlabel("r ($\AA$)")
pl.ylabel("g(r)")
if endi==None:
    endi=Ntot
pl.title("%s [%d,%d]"%(pcfile,starti,endi))
pl.show()

