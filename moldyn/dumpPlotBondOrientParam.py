#!/usr/bin/python

#mine
import plotRemote as pr #must be imported before pylab
from lammpsIO import dumpReadNext
import orderParam
from datatools import wsmooth #window smoothing
#notmine
import pylab as pl
import sys
from numpy import *

def usage():
    print "%s <lammps-dump> <l=0,1,2..> <avg=0> <cutoff=4.0> <nbins=256> <smooth=0>"%sys.argv[0]
    print "Avg=1 enables averaging accross all configurations in the dump-file"

if len(sys.argv)<2:
    usage()
    exit(0)

dumpfile=sys.argv[1]

l=int(sys.argv[2])

avg=0
if len(sys.argv)>=4:
    avg=int(sys.argv[3])

cutoff=4.
if len(sys.argv)>=5:
    cutoff=float(sys.argv[4])

nbins=256
if len(sys.argv)>=6:
    nbins=int(sys.argv[5])

smooth=0
if len(sys.argv)>=7:
    smooth=int(sys.argv[6])

dump=open(dumpfile,"r").readlines()

def plotting(rbins,rdist):
    pl.plot(rbins,rdist)
    pl.ylabel(r"$P \left (Q_%d \right )"%l)
    pl.xlabel(r"$Q_%d$"%l)
    pl.title("Bond Orientation | File:%s | Step:%s"%(sys.argv[1],head.split()[-1]))
    pr.prshow("bondOrientation.png")

if avg==1:
    cntavg=zeros(nbins)

cnt=0
while True:
    try:
        dump,bounds,types,atoms,head = dumpReadNext(dump)
        cnt+=1
    except Exception as inst:
        #last of the configurations
        if avg==0:
            print "No more configs in dump file."
        else:
            print "Found %d configurations."%cnt
        break

    bounds=[[0,bounds[0][0]],[0,bounds[1][1]],[0,bounds[2][2]]]
    Qis=orderParam.bondOrientation(array(atoms),array(bounds),l,neighbCut=cutoff,ortho=True)
    (cnts,bins) = histogram(Qis,bins=nbins,range=(0,cutoff))
    
    if avg==0:
        if smooth==1:
            smCnts = wsmooth(cnts,nbins/10)
            plotting(bins,smCnts)
        else:
            plotting(bins,cnts)
    else:
        cntavg+=cnts
if avg==1:
    plotting(bins,cntavg/cnt)
