#!/usr/bin/python

#The CHGCARs analyzed by this script need to be part of rectangular simulations
#parallepiped simulations can be changed to cubic/rectangular using the script
#poscarRectify.py
#This greatly simplifies the application of periodic boundary conditions

#This script loops over atomic bonds and pulls out the charge density between atoms
#at specific lengths (user defined).  Then plots the average density.
import sys
import operator
from scipy import *
import pylab as pl
#mine
from struct_tools import dist_periodic
import chgcarIO
from neighbors import voronoiNeighbors
from fieldPointAnalysis import fieldNeighbors1D

def chgcarBondAnalysis(chgcarfile,bondLengths,normalize=False,verbose=False):

    BLs=bondLengths
    nBLs=len(BLs)

    #Parse CHGCAR
    chgcar=open(chgcarfile,"r").readlines()
    (v1,v2,v3,atypes,axs,ays,azs,header),gridSize,chg = chgcarIO.read(chgcar)
    basis=asarray([v1,v2,v3])
    lengths=array([v1[0],v2[1],v3[2]])
    atoms=array(zip(axs,ays,azs))

    #Grid properties
    nGridPoints=reduce(operator.mul,gridSize)

    #Neighbors
    halfNeighbors=voronoiNeighbors(atoms,basis,atypes,style='half')

    a=chg.shape
    AvgChg=sum([sum([sum(line) for line in plane]) for plane in chg])/a[0]/a[1]/a[2]
    if verbose:
        print "Average CHG value:",AvgChg

    #Evaluate the CHG between each nieghbor pair
    avgchgline,xchglines,ychglines,halfNeighbors=fieldNeighbors1D(atoms,atypes,basis,chg,gridSize,halfNeighbors)

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
        avgOBL,nobl,avgIBL,nibl=chgcarBondAnalysis(chgcarfile,bondLengths,normalize,verbose)
        avgOBLs.append(avgOBL)
        nobls.append(nobl)
        avgIBLs.append(avgIBL)
        nibls.append(nibl)

    colors=["red","blue","green","purple","orange","black"]
    styles=["-","-.","--","o"]
    widths=[1,2,2,1]
    #Plotting results
    pl.figure()
    for bl in range(len(bondLengths)):
        ls=styles[bl]
        lw=widths[bl]
        #for i,(avgo,avgi) in enumerate(zip(list(avgOBLs[j]),list(avgIBLs[j]))):
        for cf in range(len(chgcarfiles)):
            c=colors[cf]
            avgi=avgIBLs[cf][bl]
            L=len(avgi)
            print "%s: <%g nm Charge at halfway point is: %g"%(chgcarfiles[cf],bondLengths[bl],(avgi[24]+avgi[25])/2.0)
            pl.plot([bondLengths[bl]*(x+0.5)/L for x in range(-L/2,L/2,1)],avgi,label="<"+str(bondLengths[bl])+" #B="+str(int(nibl[bl]))+", "+fileLabels[cf],c=c,ls=ls,lw=lw)

    pl.xlabel("Distance $\AA$")
    if normalize:
        pl.ylabel(r"Charge Density $\rho / \rho_{avg}$")
    else:
        pl.ylabel(r"Charge Density $e/\AA^3$")
    pl.legend(title="Bond Len (nm)",loc=0)
    if len(chgcarfiles)==1:
        pl.title(chgcarfile)
    else:
        pl.title(", ".join(chgcarfiles))
    pl.show()

