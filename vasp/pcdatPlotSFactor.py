#!/usr/bin/python

import pylab as pl
import sys
#mine
from structureFactor import structureFactor

pcfile=sys.argv[1]
pcdat=open(pcfile,"r").readlines()
head=pcdat[:12]
Nbins=int(pcdat[6])
delr=float(pcdat[8])*1E10

#chop off header, file is now homogeneous
pcdat=pcdat[11:]
starti=0
endi=None
xx=[i*delr for i in range(Nbins)]
yy=[map(float,pcdat[i*Nbins+i+1:(i+1)*Nbins+i+1]) for i in range(len(pcdat)/(Nbins+1))]
Ntot=len(yy)
Nused=len(yy[starti:endi])
yy1 = [sum(y)/len(y) for y in zip(*yy[starti:endi])]

aa,bb=structureFactor(xx,yy1,11.0,1024,ndens=0.045297)
pl.plot(aa,bb)
pl.show()
