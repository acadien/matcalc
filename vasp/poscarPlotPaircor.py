#!/usr/bin/python

import sys
import pylab as pl
#mine
from poscarIO import readposcar
from paircor import paircor,partpaircor
from duplicate import duplicate26
from datatools import wsmooth

def usage():
    print "%s <poscar/BestPOSCARs file> <nbins=1000> <smooth=0> <type1> <type2>"%sys.argv[0]
    print "For parital pair-correlation indicate the desired types (1,2,3...) as defined by the ordering of the POTCAR." 
    print "Note: Periodicity of the system is accounted for."

if len(sys.argv) not in [2,3,4,6]:
    usage()
    exit(0)

poscar=open(sys.argv[1],"r").readlines()

nbins=1000
if len(sys.argv)>=3:
    nbins=int(sys.argv[2])

smooth=0
if len(sys.argv)>=4:
    smooth=int(sys.argv[3])

type1=-1
type2=-1
part=0
if len(sys.argv)==6:
    type1=int(sys.argv[4])
    type2=int(sys.argv[5])
    part=1


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

    #Check for valid types
    if part==1:
        if (type1 not in types):
            print "Error: the type: %d is not in the types list.  Valid options are:%s"%(type1,thetypes)
            usage()
            exit()
        if (type2 not in types):
            print "Error: the type: %d is not in the types list.  Valid options are:%s"%(type2,thetypes)
            usage()
            exit()

    N=len(types)

    #Duplicate
    datoms,dtypes,dbasis=duplicate26(zip(ax,ay,az),types,zip(v1,v2,v3))

    #Correlate
    if part==1:
        [rbins,rdist]=partpaircor(datoms,types,type1,type2,inloop=N,nbins=nbins)
    else:
        [rbins,rdist]=paircor(datoms,inloop=N,nbins=nbins)
    
    rdist=[i/(27.0**0.5) for i in rdist]

    pl.figure()

    #Smooth
    if smooth==1:
        #rdist=windowavg(rdist,50)
        #bandpass(rbins,rdist,0.1,4.0)
        smdist=rdist[:]
        smdist = wsmooth(smdist,50)
        pl.plot(rbins,smdist)

    #Plotting
    pl.plot(rbins,rdist)
    pl.xlabel("radius (A)")
    pl.ylabel("g(r)")
    if part==1:
        pl.title("Partial Pair Correlation for types %d and %d"%(type1,type2))
    else:
        pl.title("Pair Correlation")
    pl.show()
