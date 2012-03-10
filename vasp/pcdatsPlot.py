#!/usr/bin/python

import sys
import pylab as pl
from scipy import *

def usage():
    print "Usage: %s <PCDAT1> <PCDAT2> ... <1:generate average>"%sys.argv[0]

if len(sys.argv)<2:
    usage()
    exit(0)

enableAvg=False
if sys.argv[-1]=="1":
    sys.argv.pop(-1)
    enableAvg=True

nys=list()
yys=list()
xxs=list()
pcdats=sys.argv[1:]
for pcfile in pcdats:
    pcdat=open(pcfile,"r").readlines()
    head=pcdat[:12]
    Nbins=int(pcdat[6])
    delr=float(pcdat[8])*1E10

    pcdat=pcdat[11:]

    xx=[i*delr for i in range(Nbins)]
    xxs.append(xx)

    yy=list()
    while len(pcdat)>Nbins:
        pcdat.pop(0)
        yy.append(map(float,pcdat[0:Nbins]))
        pcdat=pcdat[Nbins:]
    nys.append(len(yy))
    yys.append([sum(y)/len(y) for y in zip(*yy)])

Niter=sum(nys)
#yys=[[sum(j)/len(j) for j in zip(*i)] for i in yys]
for x,y in zip(xxs,yys):
    pl.plot(x,y)
if enableAvg:
    try:
        yavg=[sum(array(yy)*nys)/Niter for i,yy in enumerate(zip(*yys))]
        pl.plot(xxs[0],yavg,lw=4)
        a=open("PCDATAVG","w")
        a.write("".join(head))
        a.write("\n".join(map(str,yavg)))
    except NameError:
        pl.legend([str(i) for i in range(len(yy))])    
else:
    pl.legend(sys.argv[1:]+["Average"])    


pl.xlabel("r ($\AA$)")
pl.ylabel("g(r)")
if enableAvg:
    pl.title("%d Total Iterations"%(Niter))
pl.show()

