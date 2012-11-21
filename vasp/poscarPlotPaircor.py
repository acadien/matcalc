#!/usr/bin/python

import sys
import pylab as pl
from scipy import array
#mine
from poscarIO import readposcar
from paircor import paircor_periodic,partpaircor
from datatools import wsmooth

def usage():
    print "%s <poscar/BestPOSCARs file> <cutoff=10.0> <nbins=256> <smooth=0> <type1> <type2>"%sys.argv[0]
    print "For parital pair-correlation indicate the desired types (1,2,3...) as defined by the ordering of the POTCAR." 
    print "Note: Periodicity of the system is accounted for."

if len(sys.argv) not in [2,3,4,6]:
    usage()
    exit(0)

poscar=open(sys.argv[1],"r").readlines()

cutoff=10.
if len(sys.argv)>=3:
    cutoff=float(sys.argv[2])

nbins=256
if len(sys.argv)>=4:
    nbins=int(sys.argv[3])

smooth=0
if len(sys.argv)>=5:
    smooth=int(sys.argv[4])

type1=-1
type2=-1
part=0
if len(sys.argv)==7:
    type1=int(sys.argv[5])
    type2=int(sys.argv[6])
    part=1


while True:
    [v1,v2,v3,atypes,ax,ay,az,head,poscar] = readposcar(poscar)

    atoms=array(zip(ax,ay,az))

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
    basis=array([v1,v2,v3])
    #Duplicate
    #datoms,dtypes,dbasis=duplicate26(zip(ax,ay,az),types,zip(v1,v2,v3))

    #Correlate
    if part==1:
        [rbins,rdist]=partpaircor(atoms,types,type1,type2,inloop=N,nbins=nbins)
    else:
        print atoms
        [rbins,rdist]=paircor_periodic(atoms,basis,cutoff=cutoff,nbins=nbins)
    
    rdist=[i for i in rdist]

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
