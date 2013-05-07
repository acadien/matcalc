#!/usr/bin/python

#mine
import plotRemote as pr
import poscarIO
from duplicate import duplicate26
from datatools import wsmooth
from struct_tools import dist,dist_periodic
from neighbors import neighbors
from rdf import adf
#notmine
import sys
import pylab as pl
from scipy import array
from mpl_toolkits.mplot3d import Axes3D
from math import fabs
from numpy import array,zeros
def usage():
    print "%s <r-cut> <poscar files>"%sys.argv[0]
    print "Note: Periodicity of the system is accounted for."

if len(sys.argv) < 3:
    usage()
    exit(0)

poscars=[open(sys.argv[i],"r").readlines() for i in range(2,len(sys.argv))]

nbins=360
rcut=float(sys.argv[1])
rdists=zeros(nbins)
for poscar in poscars:
    [basis,atypes,atoms,head,poscar] = poscarIO.read(poscar)
    atoms=array(atoms)

    N=atoms.shape[0]

    #Neighbors needed for angular distribution
    lengths=array([basis[0][0],basis[1][1],basis[2][2]])
    bounds = [[0,lengths[0]],[0,lengths[1]],[0,lengths[2]]]
    neighbs = neighbors(atoms,bounds,rcut)
    [rbins,rdist]=adf(atoms,neighbs,basis,nbins=nbins)
    rdists+=rdist
rdists/=len(poscars)
#Smooth
smdist=rdists[:]
smdist = wsmooth(smdist,10)

    #Plotting
pl.figure()
pl.plot(rbins,rdists)
pl.plot(rbins,smdist)

pl.xlabel("Angle (deg)")
pl.ylabel("Count")
pl.title("Distribution of Angles %s"%sys.argv[1])
pr.prshow("poscarADF.png")
    
