#!/usr/bin/python

#mine
import plotRemote as pr #must be imported before pylab
from lammpsIO import dumpReadNext
from rdf import rdf_periodic
from datatools import wsmooth #window smoothing
#notmine
import pylab as pl
import sys
from numpy import array,zeros

def usage():
    print "%s <lammps-dump> <avg=0> <cutoff=10.0> <nbins=256> <smooth=0>"%sys.argv[0]
    print "Avg>=1 enables averaging, using at most the number of averages requested."

if len(sys.argv)<2:
    usage()
    exit(0)

dumpfile=sys.argv[1]

avg=0
if len(sys.argv)>=3:
    avg=int(sys.argv[2])

cutoff=10.
if len(sys.argv)>=4:
    cutoff=float(sys.argv[3])

nbins=256
if len(sys.argv)>=5:
    nbins=int(sys.argv[4])

smooth=0
if len(sys.argv)>=6:
    smooth=int(sys.argv[5])

dump=open(dumpfile,"r").readlines()

def plotting(rbins,rdist):
    pl.plot(rbins,rdist)
    pl.xlabel("radius ($\AA$)")
    pl.ylabel("g(r)")
    pl.title("RDF | File:%s | Step:%s"%(sys.argv[1],head.split()[-1]))
    pr.prshow("RDF_lammps.png")

if avg>=1:
    rdistavg=zeros(nbins)

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
    
    [rbins,rdist]=rdf_periodic(array(atoms),array(bounds),cutoff=cutoff,nbins=nbins)
    
    if avg==0:
        if smooth==1:
            smdist=rdist[:]
            smdist = wsmooth(smdist,nbins/10)
            plotting(rbins,smdist)
        else:
            plotting(rbins,rdist)
    else:
        rdistavg+=rdist

    if avg>=1 and cnt%avg==0 and cnt>0: 
        plotting(rbins,rdistavg/avg)
        rdistavg=zeros(nbins)


