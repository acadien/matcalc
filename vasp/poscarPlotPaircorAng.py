#!/usr/bin/python

#mine
import plotRemote as pr
import poscarIO
from duplicate import duplicate26
from datatools import wsmooth
from struct_tools import neighbors,dist,dist_periodic
from voronoiNeighbors import voronoiNeighbors
from rdf import adf
#notmine
import sys
import pylab as pl
from scipy import array
from mpl_toolkits.mplot3d import Axes3D
from math import fabs

def usage():
    print "%s <poscar/BestPOSCARs file> <nbins=360> <optional:minBondLen,maxBondLen>"%sys.argv[0]
    print "Note: Periodicity of the system is accounted for."

if len(sys.argv) not in [2,3,4]:
    usage()
    exit(0)

poscar=open(sys.argv[1],"r").readlines()

nbins=360
if len(sys.argv)>=3:
    nbins=int(sys.argv[2])
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

while True:
    [v1,v2,v3,atypes,ax,ay,az,head,poscar] = poscarIO.read(poscar)

    if v1==v2==v3==-1:
        break

    j=1
    types=list()
    for i in atypes:
        types+=[j]*i
        j+=1
    thetypes=str(set(types))

    N=len(types)
    atoms=array(zip(ax,ay,az))
    basis=array([v1,v2,v3])

    #Duplicate
    #datoms,dtypes,dbasis=duplicate26(atoms,types,basis)

    #Neighbors needed for angular distribution

    lengths=array([basis[0][0],basis[1][1],basis[2][2]])

    #dbounds=[[0.0,dbasis[0][0]],[0.0,dbasis[1][1]],[0.,dbasis[2][2]]]
    #neighbs=voronoiNeighbors(atoms=datoms,atypes=dtypes,basis=dbasis,style='full')
    neighbs=voronoiNeighbors(atoms=atoms,atypes=atypes,basis=basis,style='full')

    if bl!=-1 and bw!=-1:
        neighbs=[[j for j in neighbs[i] if fabs(dist_periodic(atoms[i],atoms[j],lengths)-bl)<bw] for i in range(N)]
    #Correlate
    [rbins,rdist]=adf(atoms,neighbs,basis,nbins=nbins)
    
    #Smooth
    #if smooth==1:
    #    smdist=rdist[:]
    #    smdist = wsmooth(smdist,50)

    #Plotting
    pl.figure()
    #if smooth==1:
    #    pl.plot(rbins,smdist)
    #    pl.plot(rbins,rdist)
    #else:
    pl.plot(rbins,rdist)

    pl.xlabel("Angle (deg)")
    pl.ylabel("Count")
    pl.title("Distribution of Angles %s"%sys.argv[1])
    pr.prshow("poscarADF.png")
    
