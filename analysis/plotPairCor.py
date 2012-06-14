#!/usr/bin/python

import sys
import pylab as pl
#mine
from paircorIO import readPairCor

def usage():
    print "Usage: %s <pair correlation files, space seperated> <N=normalize>"%sys.argv[0].split("/")[-1]

norm=False
if sys.argv[-1] in ["n","N","normalize","Normalize"]:
    norm=True
    fnames=sys.argv[1:-1]
else:
    fnames=sys.argv[1:]

for fname in fnames:
    header,bins,vals=readPairCor(fname)
    
    if norm:
        tvals=sum(vals)
        vals=[i/tvals for i in vals]
        pl.yticks([])
    
    pl.plot(bins,vals,label=" ".join(fname.strip(".data").split("_")[1:]))

pl.xlabel("r $(\AA)$")
pl.ylabel("Bond Count")
if len(fnames)>1:
    pl.legend(loc=0)
else:
    pl.title("%s"%header)
pl.show()
