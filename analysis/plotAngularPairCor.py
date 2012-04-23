#!/usr/bin/python

import sys
import pylab as pl
#mine
from angularpaircorIO import readAngularPairCor

def usage():
    print "Usage: %s <angular pair correlation files, space seperated> <N=normalize>"

norm=False
if sys.argv[-1] in ["n","N","normalize","Normalize"]:
    norm=True
    fnames=sys.argv[1:-1]
else:
    fnames=sys.argv[1:]

for fname in fnames:
    header,minlen,maxlen,bins,vals=readAngularPairCor(fname)
    if norm:
        tvals=sum(vals)
        vals=[i/tvals for i in vals]
        pl.yticks([])
    pl.plot(bins,vals,label=" ".join(fname.strip(".data").split("_")[1:]))

pl.xlabel("Angle")
pl.ylabel("Number Bonds Triples")
if len(fnames)>1:
    pl.legend(loc=0)
else:
    pl.title("%s %g to %g"%(header,minlen,maxlen))
pl.show()
