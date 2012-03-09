#!/usr/bin/python

import sys
import pylab as pl
#mine
from poscarIO import readposcar
from duplicate import duplicate26
from datatools import wsmooth
from struct_tools import neighbors
from paircor import paircor_ang

def usage():
    print "%s <poscar/BestPOSCARs file> <nbins=1000> <smooth=0>"%sys.argv[0]
    print "Note: Periodicity of the system is accounted for."

if len(sys.argv) not in [2,3,4]:
    usage()
    exit(0)

poscar=open(sys.argv[1],"r").readlines()

nbins=1000
if len(sys.argv)>=3:
    nbins=int(sys.argv[2])

smooth=0
if len(sys.argv)>=4:
    smooth=int(sys.argv[3])

while True:
    [v1,v2,v3,atypes,ax,ay,az,head,poscar] = readposcar(poscar)

    if v1==v2==v3==-1:
        break

    j=1
    types=list()
    for i in atypes:
        types+=[j]*i
        j+=1
    thetypes=str(set(types))

    N=len(types)
    atoms=zip(ax,ay,az)

    #Duplicate
    datoms,dtypes,dbasis=duplicate26(atoms,types,zip(v1,v2,v3))

    #Neighbors needed for angular distribution
    neighbs=neighbors(datoms,5.0,style="full")

    #Correlate
    [rbins,rdist]=paircor_ang(datoms,neighbs,nbins=nbins,inloop=N)
    
    #rdist=[i/(27.0**0.5) for i in rdist]
    #Smooth
    if smooth==1:
        #rdist=windowavg(rdist,50)
        #bandpass(rbins,rdist,0.1,4.0)
        smdist=rdist[:]
        smdist = wsmooth(smdist,50)

    #Plotting
    pl.figure()
    if smooth==1:
        pl.plot(rbins,smdist)
        pl.plot(rbins,rdist)
    else:
        pl.plot(rbins,rdist)

    pl.xlabel("Angle (deg)")
    pl.ylabel("Count")
    pl.title("Distribution of Angles")
    pl.show()
