#!/usr/bin/python

import sys
import pylab as pl
#mine
from angularpaircorIO import readAngularPairCor

def usage():
    print "Usage: %s <angular pair correlation file>"

fname=sys.argv[1]
header,minlen,maxlen,bins,vals=readAngularPairCor(fname)
pl.plot(bins,vals)
pl.title("%s %g to %g"%(header,minlen,maxlen))
pl.xlabel("Angle")
pl.ylabel("Number Bonds Triples")
pl.show()
