#!/usr/bin/python

#CHGCAR is read in, along with a desired bond length to look for
#

#This script loops over atomic bonds and pulls out the charge density between atoms
#at specific lengths (user defined).  Then plots the average density.
import sys
import operator
from scipy import *
import pylab as pl
#mine
from struct_tools import dist_periodic
from chgcarIO import readchgcar
from voronoiNeighbors import voronoiNeighbors
from fieldPointAnalysis import fieldNeighbors3D

def chgcarBondAnalysis(chgcarfile,bondLengths,normalize=False,verbose=False):

    BLs=bondLengths
    nBLs=len(BLs)

    #Parse CHGCAR
    chgcar=open(chgcarfile,"r").readlines()
    (v1,v2,v3,atypes,axs,ays,azs,header),gridSize,chg = readchgcar(chgcar)
    basis=asarray([v1,v2,v3])
    lengths=array([v1[0],v2[1],v3[2]])
    atoms=array(zip(axs,ays,azs))

    #Grid properties
    nGridPoints=reduce(operator.mul,gridSize)

    #Neighbors
    halfNeighbors=voronoiNeighbors(atoms=atoms,basis=basis,atypes=atypes,style='half')

    a=chg.shape
    AvgChg=sum([sum([sum(line) for line in plane]) for plane in chg])/a[0]/a[1]/a[2]
    if verbose:
        print "Average CHG value:",AvgChg

    #Evaluate the CHG between each nieghbor pair
    avgchgline,xchglines,ychglines,halfNeighbors=fieldNeighbors3D(atoms,atypes,basis,chg,gridSize,halfNeighbors)

    #Cutoff neighbors that fall below the thresh-hold
    Ninterp=len(avgchgline)
    avgIBL=zeros([nBLs,Ninterp]) #avg line inside the bond length
    avgOBL=zeros([nBLs,Ninterp]) #avg line outside the bond length
    nibl=zeros(nBLs)
    nobl=zeros(nBLs)

    cnt=0
    for a,jNeighbors in enumerate(halfNeighbors):
        atoma=atoms[a]
        for b,jNeighb in enumerate(jNeighbors):
            atomb=atoms[jNeighb]
            d=dist_periodic(atoma,atomb,lengths)
            vals=ychglines[a][b]
            xx=xchglines[a][b]
            for j,bl in enumerate(BLs):
                if normalize:
                    avals=array(vals)/AvgChg
                else:
                    avals=array(vals)
                if d<bl:
                    nibl[j]+=2
                    avgIBL[j]+=avals
                    avgIBL[j]+=avals[::-1]
                else:
                    cnt+=1
                    nobl[j]+=2
                    avgOBL[j]+=avals
                    avgOBL[j]+=avals[::-1]
    for i in range(nBLs):
        if nibl[i]==0:
            nibl[i]=1
        if nobl[i]==0:
            nobl[i]=1
    avgOBL=[avgOBL[i]/nobl[i] for i in range(nBLs)]
    avgIBL=[avgIBL[i]/nibl[i] for i in range(nBLs)]
    return avgOBL,nobl,avgIBL,nibl

if __name__=="__main__":
    def usage():
        print "Usage:"
        print "%s <CHGCAR0,CHGCAR1...> <maxL0,maxL1,maxL2...> <[N,n]ormalize by total avg charge>"%(sys.argv[0])
        print "Note, CHGCAR must be part of a cubic/rectangular simulation"

    if len(sys.argv)<3:
        usage()
        exit(0)

    chgcarfiles=sys.argv[1].split(",")
    bondLengths=map(float,sys.argv[2].split(","))
    normalize=False
    if len(sys.argv)>3:
        if sys.argv[3][0] in ["n","N"]:
            normalize=True
    verbose=True


    try:
        fileLabels=[i.split("_")[1] for i in chgcarfiles]
    except:
        fileLabels=chgcarfiles

    avgOBLs=list()
    nobls=list()
    avgIBLs=list()
    nibls=list()
    for chgcarfile in chgcarfiles:
        avgGrid,neighbs=chgcarBondAnalysis(chgcarfile,bondLengths,normalize,verbose)
        #WRITE THIS TO FILE
        print avgGrid


