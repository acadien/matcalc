#!/usr/bin/python

import sys
import pylab as pl
#mine
from angularpaircorIO import readAngularPairCor

def usage():
    print "Usage: %s <angular pair correlation files, space seperated>"

fnames=sys.argv[1:]
for fname in fnames:
    header,minlen,maxlen,bins,vals=readAngularPairCor(fname)
    pl.plot(bins,vals,label=" ".join(fname.strip(".data").split("_")[1:]))

pl.xlabel("Angle")
pl.ylabel("Number Bonds Triples")
if len(fnames)>1:
    pl.legend(loc=0)
else:
    pl.title("%s %g to %g"%(header,minlen,maxlen))
pl.show()
