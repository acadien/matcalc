#!/usr/bin/python
import sys
from math import *
from scipy import array
import pylab as pl
#mine
from paircor import paircor_ang
from voronoiNeighbors import *

#Calculates the pair-correlation function for a set of atoms in an XDATCAR file and writes it to a file

def usage():
    print "Usage: %s <Outcar> <optional:numbins=360> <optional:minBondLen,maxBondLen>"%sys.argv[0]
    print "If a bondlen is given only bonds of length bondeln+/-err will be evaluated"

def outcarPairCorAng(outcarfile,nbins,bl=-1.0,bw=-1.0):
    outcar=open(outcarfile,"r")
    tbinvals=list()

    #Grab ion types
    while True:
        line=outcar.readline()
        if "ions per type" in line:
            break
    atypes=map(int,line.split("=")[1].split())

    #Grab basis vectors
    while True:
        line=outcar.readline()
        if "direct lattice vectors" in line:
            break
    basis=array([map(float,outcar.readline().split()[:3]) for i in range(3)])
    lengths=array([basis[0][0],basis[1][1],basis[2][2]])

    #Grab atom positions and perform paircorang analysis
    count=0
    posit=False
    for line in outcar.readlines():
        if posit:
            if "--" in line:
                if len(atoms)==0:
                    continue
                else:
                    #Analysis
                    atoms=array(atoms)
                    neighbs=voronoiNeighbors(atoms=atoms,basis=basis,atypes=atypes,style='full')
                    if bl!=bw!=-1:
                        specbonds=[list() for i in range(len(atoms))]
                        for i in range(len(atoms)):
                            ai=dot(atoms[i],basis)
                            for j in neighbs[i]:
                                aj=dot(atoms[j],basis)
                                d=dist(ai,aj)
                                if d!=dist_periodic(ai,aj,lengths):
                                    for c in range(3):
                                        d = aj[c]-ai[c]
                                        if d>lengths[c]/2.0: d -= lengths[c]
                                        elif d<-lengths[c]/2.0: d += lengths[c]
                                if fabs(d-bl)<bw:
                                    specbonds[i].append(j)
                    else:
                        specbonds=neighbs
                    [angs,abins]=paircor_ang(atoms,specbonds,basis,nbins=nbins)
                    tbinvals.append(abins)
                    print count
                    if count > 1000:
                        #tbinvals/=float(count)
                        break
                    posit=False
            else:
                atoms.append(map(float,line.split()[:3]))
        elif "POSITION" in line:
            atoms=list()
            posit=True
            count+=1
    
    return angs,tbinvals

if __name__=="__main__":
    if len(sys.argv)<2:
        usage()
        exit(0)
    bl=-1.
    bw=-1.
    nbins=360
    if len(sys.argv)>=3:
        nbins=int(sys.argv[2])
        if len(sys.argv)==4:
            bmin,bmax=map(float,sys.argv[3].split(","))
        bl=(bmin+bmax)/2.
        bw=bmax-bl

    outcarfile=sys.argv[1]

    angs,tbinvals=outcarPairCorAng(outcarfile,nbins,bl,bw)
    pl.plot(angs,map(lambda x:sum(x)/len(x),zip(*tbinvals)))
    if bl==bw==-1.:
        pl.title("Angular Distribution of %s"%sys.argv[1])
    else:
        pl.title("Angular Distribution of %s %g$\pm$%g"%(sys.argv[1],bl,bw))
    pl.show()


