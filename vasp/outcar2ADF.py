#!/usr/bin/python
import sys
from math import *
from scipy import array
import pylab as pl
#mine
from rdf import adf
from neighbors import voronoiNeighbors
import adfIO

def usage():
    print "Usage: %s <Outcar> <optional:numbins=360> <optional:minBondLen,maxBondLen>"%sys.argv[0]
    print "If a bondlen is given only bonds with a length between minlen and maxlen will be evaluated."

def outcarADF(outcarfile,nbins,bl=-1.0,bw=-1.0):
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

    #Grab atom positions and perform adf calculation
    count=0
    posit=False
    for line in outcar:
        if posit:
            if "--" in line:
                if len(atoms)==0:
                    continue
                else:
                    #Analysis
                    if max(zip(*atoms)[0])<1.01:
                        atoms=dot(array(atoms),basis)
                    else:
                        atoms=array(atoms)
                    neighbs=voronoiNeighbors(atoms,basis,atypes,style='full')
                    if bl!=-1 and bw!=-1:
                        specbonds=[list() for i in range(len(atoms))]
                        for i in range(len(atoms)):
                            ai=atoms[i]#dot(atoms[i],basis)
                            for j in neighbs[i]:
                                aj=atoms[j]#dot(atoms[j],basis)
                                d=dist_periodic(ai,aj,lengths)
                       
                                if fabs(d-bl)<bw:
                                    if d<1.0:
                                        print d
                                        print dist_periodic(ai,aj,lengths)
                                        print ai
                                        print aj
                                        exit(0)
                                    specbonds[i].append(j)
                    else:
                        specbonds=neighbs
                    [angs,abins]=adf(atoms,specbonds,basis,nbins=nbins)
                    tbinvals.append(abins)
                    print count
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
    bmin=0.0
    bmax=10.0
    if len(sys.argv)>=3:
        nbins=int(sys.argv[2])
    if len(sys.argv)>=4:
        if len(sys.argv)==4:
            bmin,bmax=map(float,sys.argv[3].split(","))
    bl=float(bmin+bmax)/2.
    bw=bmax-bl

    outcarfile=sys.argv[1]
    
    try:
        num=sys.argv[1].split("_")[1]
    except IndexError:
        num=None

    angs,tbinvals=outcarADF(outcarfile,nbins,bl,bw)
    avgvals=map(lambda x:sum(x)/len(x),zip(*tbinvals))

    if num:
        adffile="ADF_%s_%g-%g.data"%(num,bmin,bmax)
    else:
        adffile="ADF_%g-%g.data"%(bmin,bmax)

    header="Angular Distribution from file %s."%(sys.argv[0])

    print "Writing %s."%apcfile
    adfIO.write(adffile,header,bmin,bmax,angs,avgvals)
