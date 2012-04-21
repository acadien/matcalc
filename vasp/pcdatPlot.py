#!/usr/bin/python

import sys
import pylab as pl
from scipy import *
#mine
from datatools import gaussSmooth

def usage():
    print "Usage: %s <PCDAT file> <optional:start iter(0)> <optional:end iter(-1)> <1:write average>"%sys.argv[0].split("/")[-1]

starti=0
endi=None
enableAvg=False
if len(sys.argv)<2:
    usage()
    exit(0)
elif len(sys.argv)>3:
    starti=int(sys.argv[2])
    endi=int(sys.argv[3])
if len(sys.argv) in [3,5]:
    if sys.argv[2]=='1' or sys.argv[4]=='1':
        enableAvg=True

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
yy1 = [sum(y)/len(y) for y in zip(*yy[starti:endi])]
#yys = gaussSmooth([sum(y)/len(y) for y in zip(*yy[starti:endi])],20,mode='valid')

print "Using %d iterations of %d possible"%(Nused,Ntot)

if enableAvg:
    if endi==None:
        a=open("PCDATAVG_%d_%d"%(0,Nused),"w")
    else:
        a=open("PCDATAVG_%d_%d"%(starti,endi),"w")
    a.write("".join(head))
    a.write("\n".join(map(str,yy1)))
    a.close()

pl.plot(xx,yy1)
pl.xlabel("r ($\AA$)")
pl.ylabel("g(r)")
if endi==None:
    endi=Ntot
pl.title("%s [%d,%d]"%(pcfile,starti,endi))
pl.show()

