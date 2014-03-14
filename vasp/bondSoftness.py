#!/usr/bin/python

import plotRemote as pr
import sys
import pylab as pl
import numpy as np
from neighbors import neighbors
#mine
import procarIO
import poscarIO
from struct_tools import dist_periodic

#Evaluates the fukui softness functions using the local density of states

def eFermiParse(doscarFile):
    dcarF=open(doscarFile,"r")
    for i in range(5): dcarF.next() 
    return float(dcarF.next().split()[3])

def usage():
    print "procarPlot.py <PROCAR File> <DOSCAR File> <CONTCAR File> <Neighbors File>"
    print "Neighbors file is best generated by elfcarNeighbors.py"
    print "To generate the correct PROCAR make sure to have the following settings in your INCAR:\nISMEAR=0(semiconductors), 1(Metals+SIGMA)\nPREC=high\nLORBIT=11\nNPAR=1"

if len(sys.argv)!=5:
    usage()
    exit(0)

procarFile=sys.argv[1]
doscarFile=sys.argv[2]
poscarFile=sys.argv[3]
neighbsFile=sys.argv[4]

#Fermi Energy From DOSCAR
eFermi = eFermiParse(doscarFile)

#Band's and Band Energy from Procar
nIon,nGridPoints,bandGrid,occGrid = procarIO.readLDOS(procarFile)
bandGrid-=eFermi

#Ion positions from CONTCAR/POSCAR
basis,atypes,atoms,head,poscar = poscarIO.read(poscarFile)
atoms=np.asarray(map(np.asarray,atoms))
lengths=np.asarray(map(np.linalg.norm,basis))
bounds=[[0,lengths[0]],[0,lengths[1]],[0,lengths[2]]]

#Parse Neighbors file
try:
    rcut=float(neighbsFile)
    neighbs=neighbors(atoms,bounds,rcut,style="full")
except ValueError:
    neighbs=[map(int,i.split()) for i in open(neighbsFile,"r").readlines()[1:]]

#chemical potential integral
delU = 0.1
uBounds = np.where(np.logical_and(bandGrid > -delU, bandGrid < 0))[0]
bandBounds=bandGrid[uBounds]
d= (bandBounds.max()-bandBounds.min())/len(bandBounds)

rlen=list()
soft=list()
for i in range(nIon):
    for j in neighbs[i]:
        r=dist_periodic(atoms[i],atoms[j],lengths)
        
        #g(E,atomi)
        ga=occGrid[i][uBounds]
        gb=occGrid[j][uBounds]
        #soft.append(sum(ga*gb)/d/delU*2)
        soft.append(sum(ga+gb)/delU)
        rlen.append(r)

pl.scatter(rlen,soft)
pl.show()