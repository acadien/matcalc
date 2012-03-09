#!/usr/bin/python
import sys,os

if __name__=="__main__":
    run_remote=False
    if len(sys.argv)>=4:
        if sys.argv[3]=='1':
            run_remote=True
            import matplotlib
            matplotlib.use('Agg') 
            figname=sys.argv[4]

import operator
from scipy import *
#mine
from elfcarAnalyze import elfNeighbors
from voronoiNeighbors import voronoiNeighbors
from struct_tools import *
from poscarIO import *

if __name__=="__main__":

    def usage():
        print "Usage:"
        print "%s <ELFCAR> <POSCAR/\'None\'> <1:remote> <figure name>"%(sys.argv[0])
        print "Note POSCAR & ELFCAR must be part of a cubic/rectangular simulation"
        print "If POSCAR is given, the voronoi neighbors are generated for comparison with ELFCAR."

    if len(sys.argv)<2:
        usage()
        exit(0)

    #Parse ELFCAR
    elfcar=open(sys.argv[1],"r").readlines()
    v1,v2,v3,atypes,axs,ays,azs,header,elfcar = readposcar(elfcar)
    basis=asarray([v1,v2,v3])
    bounds=[[0.,v1[0]],[0.,v2[1]],[0.,v3[2]]]
    atoms=asarray(zip(axs,ays,azs))
    nAtoms=len(atoms)
    r=5.0 #5 Angstroms by default

    #Generate Voronoi Neighbor list
    voronoi_enable=False
    if sys.argv[2] not in ["none","n","N","None","NONE"]:
        poscar=open(sys.argv[2],"r").readlines()
        voronoi_enable=True
        vFullNeighbors=voronoiNeighbors(poscar=poscar)
        vHalfNeighbors=voronoiNeighbors(poscar=poscar,style='half')

    #Read in Grid properties
    elfcar.pop(0)
    gridSize=[int(i) for i in elfcar.pop(0).split()]
    nGridPoints=reduce(operator.mul,gridSize)

    #Read in ELF data
    elf=asarray([float(i) for i in ("".join(elfcar)).split()])
    elf.shape=gridSize

    Nneighbs=30 #max number of neighbors

    """
    #Get the (raw) processed neighbor list
    rpHalfNeighbors,avgelf,avgfail=elfNeighbors(atoms,basis,gridSize,elf,r,bondCutoff=0.3)    

    #Raw neighbor list
    halfNeighbors=neighbors(atoms,bounds,r,style='half')
    fullNeighbors=map(list,halfNeighbors)
    for iNeighb,neighbs in enumerate(halfNeighbors):
        for jNeighb in neighbs:
            fullNeighbors[jNeighb].append(iNeighb)
    
#    print "Number of Neighbors before elimination:",sum([len(i) for i in fullNeighbors])
    #=====================
    #The raw neighbor list
    neighbhist=zeros(Nneighbs)
    for iNeighb,neighbs in enumerate(fullNeighbors):
        ncnt=len(neighbs)
        neighbhist[ncnt]+=1
    pl.plot(neighbhist)

    rpFullNeighbors=map(list,rpHalfNeighbors)
    for iNeighb,neighbs in enumerate(rpHalfNeighbors):
        for jNeighb in neighbs:
            rpFullNeighbors[jNeighb].append(iNeighb)

#   print "Number of Neighbors after elimination:",sum([len(i) for i in fullNeighbors])
    neighbhist=zeros(Nneighbs)
    for iNeighb,neighbs in enumerate(rpFullNeighbors):
        ncnt=len(neighbs)
        neighbhist[ncnt]+=1
    pl.plot(neighbhist)
    """

    #=====================
    #Voronoi neighbor list
    vneighbhist=zeros(Nneighbs)
    for iNeighb,neighbs in enumerate(vFullNeighbors):
        ncnt=len(neighbs)
        vneighbhist[ncnt]+=1
    pl.plot(vneighbhist)

    #Get the (voronoi) processed neighbor list
    vpHalfNeighbors,ae,af=elfNeighbors(atoms,basis,atypes,gridSize,elf,halfNeighbors=vHalfNeighbors,bondCutoffs=0.3)   
    vpFullNeighbors=map(list,vpHalfNeighbors)
    for iNeighb,neighbs in enumerate(vpHalfNeighbors):
        for jNeighb in neighbs:
            vpFullNeighbors[jNeighb].append(iNeighb)
 
    vpneighbhist=zeros(Nneighbs)
    for iNeighb,neighbs in enumerate(vpFullNeighbors):
        ncnt=len(neighbs)
        vpneighbhist[ncnt]+=1
    pl.plot(vpneighbhist)
    
    #Plotting
    pl.legend(["Voronoi","Proc Voron"])
    pl.xlabel("# of Neighbors")
    pl.ylabel("Count")
    pl.title("/".join(os.path.abspath(os.curdir).split("/")[-2:]) + " Neighbor List Histogram")
#    pl.xlim([0,25])
    if run_remote:
        print figname
        pl.savefig(figname)
    else:
        pl.show()
    
