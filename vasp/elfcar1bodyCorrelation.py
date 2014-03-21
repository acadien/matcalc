#!/usr/bin/python

import chgcarIO,elfcarIO
import sys
import pylab as pl
import numpy as np
def usage():
    print "%s CHGCARfiles"%sys.argv[0].split("/")[-1]

if len(sys.argv)<2:
    usage()
    exit(0)

chgelfcarfiles = sorted(sys.argv[1:],key=lambda x:str(x.split("_")[-1]))
m=0
for chgelfcarfile in chgelfcarfiles:
    car=open(chgelfcarfile,"r").readlines()
    if "CHG" in chgelfcarfile:
        pcar,field = chgcarIO.read(car)
    if "ELF" in chgelfcarfile:
        pcar,field = elfcarIO.read(car)
    (basis,atypes,atoms,head)=pcar
    volume=np.dot(np.cross(basis[0],basis[1]),basis[2])/len(atoms)
    field=field.ravel()
    field=map(float,field.tolist())
    #field.sort()
    vals,bins,dummy=pl.hist(field,100,alpha=0.5,visible=False)
    s=sum(vals)
    vals=[v/s*100 for v in vals]
    m=max(m,max(vals))
    bins=map(lambda x:(x[0]+x[1])/2.,zip(bins[:-1],bins[1:]))
    pl.plot(bins,vals,label="Step: %s"%chgelfcarfile.split("_")[-1].split("/")[0]+"   Vol:%4.4f"%volume)
#pl.xticks([x/20.*2.5 - 0.5 for x in range(20)])
pl.ylabel("Percent Volume")
pl.xlabel("Charge (meV)")
pl.ylim([0,m*1.1])
pl.legend()
pl.show()
