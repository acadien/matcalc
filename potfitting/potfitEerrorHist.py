#!/usr/bin/python

import sys
import pylab as pl

def usage():
    print "%s <energy file>"%sys.argv[0].split("/")[-1]

if len(sys.argv)<2:
    usage()
    exit(0)

efile=open(sys.argv[1],"r").readlines()[2:]

absdelE=map(lambda x:float(x.split()[5]),efile)
pl.hist(absdelE,20)
pl.xlabel("$|\Delta E|$")
pl.ylabel("count")


print "Worst fits: \ncnt config#    |delE|"
cnt=0
mxerr=0.1
for i,val in enumerate(reversed(sorted(absdelE))):
    if val<mxerr:
        break
    print i," ",absdelE.index(val)," \t",val
    cnt=i
print "%d configurations with error above %f"%(cnt,mxerr)

pl.savefig("/home/acadien/Dropbox/GeData/hist_energy_error.png")
