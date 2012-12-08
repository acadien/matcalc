#!/usr/bin/python

#CHGCAR is read in, along with a desired bond length to look for
#a plane surrounding the atomic pairs is extracted and averaged over.

#This script loops over atomic bonds and pulls out the charge density between atoms
#at specific lengths (user defined).  Then plots the average density.
import sys
import operator
from scipy import *
import pylab as pl
#mine
from struct_tools import dist_periodic
import chgcarIO
from voronoiNeighbors import voronoiNeighbors
from fieldPointAnalysis import fieldNeighbors3D

def chgcarBondAnalysis(chgcarfile,bondLengths,normalize=False,verbose=False,Ninterps=None,loc=None):

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
    halfNeighbors=voronoiNeighbors(atoms=atoms,basis=basis,atypes=atypes,style='half')

    a=chg.shape
    AvgChg=sum([sum([sum(line) for line in plane]) for plane in chg])/a[0]/a[1]/a[2]
    if verbose:
        print "Average CHG value:",AvgChg

    #Evaluate the CHG between each nieghbor pair
    if Ninterps==None:
        Ninterps=array([1,50,50])
        avgGrids,grids,bondcounts=fieldNeighbors3D(atoms,atypes,basis,chg,gridSize,halfNeighbors,cutoffs=bondLengths,loc=loc)
    else:
        Ninterps=array([1]+Ninterps)
        avgGrids,grids,bondcounts=fieldNeighbors3D(atoms,atypes,basis,chg,gridSize,halfNeighbors,cutoffs=bondLengths,Ninterps=Ninterps,loc=loc)

    #Normalize if necessary
    if normalize:
        avgGrids=[avgGrid/AvgChg for avgGrid in avgGrids]
            
    return avgGrids,grids,bondcounts

    """
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
    """

if __name__=="__main__":
    def usage():
        print "Usage:"
        print "%s <CHGCAR0,CHGCAR1...> <maxL0,maxL1,maxL2...> <[N,n]ormalize by total avg charge> <Ninterps: x,y>"%(sys.argv[0])
        print "Note, CHGCAR must be part of a cubic/rectangular simulation"
        print "Note that the y-axis is along the bond in question"

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

    Ninterps=[50,50]
    if len(sys.argv)>4:
        Ninterps=map(int,sys.argv[5].split(","))

    for chgcarfile in chgcarfiles:

        avgGrids,grids,bondCounts=chgcarBondAnalysis(chgcarfile,bondLengths,normalize,verbose,Ninterps=Ninterps,loc="center")
        chgcar=open(chgcarfile,"r").readlines()
        poscardata,gridSize,chg = chgcarIO.read(chgcar)
        
        """
        for avgGrid,bondCount,cutoff in zip(avgGrids,bondCounts,bondLengths):
            #vtkfname=chgcarfile+"_cut%2.2f_NBond%d.vtk"%(cutoff,bondCount)
            #chgcar.writeVTK(vtkfname,poscardata,avgGrid.shape,avgGrid)
            pl.imshow(avgGrid[0])
            pl.show()
        """
        for i,grid in enumerate(grids):
            pl.savetxt("grid%d_%s.data"%(i,chgcarfile),grid)
            pl.imshow(grid)
            pl.show()
