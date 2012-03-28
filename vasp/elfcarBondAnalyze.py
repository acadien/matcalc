#!/usr/bin/python

#The ELFCAR's analyzed by this script need to be part of rectangular simulations
#parallepiped simulations can be changed to cubic/rectangular using the script
#poscarRectify.py
#This greatly simplifies the application of periodic boundary conditions

import sys
import operator
from scipy import *
import pylab as pl
#mine
from elfcarIO import readelfcar
from voronoiNeighbors import voronoiNeighbors
from fieldPointAnalysis import fieldNeighbors


def usage():
    print "Usage:"
    print "%s <ELFCAR> <bond-ELF-cutoff=0.5>"%(sys.argv[0])
    print "Note ELFCAR must be part of a cubic/rectangular simulation"

if len(sys.argv)<2:
    usage()
    exit(0)

bondCutoffs=[float(i)/10 for i in range(2,7)]
if len(sys.argv)==3:
    bondCutoffs=[float(sys.argv[2])]

#Parse ELFCAR
elfcar=open(sys.argv[1],"r").readlines()
(v1,v2,v3,atypes,axs,ays,azs,header),gridSize,elf = readelfcar(elfcar)
basis=asarray([v1,v2,v3])
bounds=[[0.,v1[0]],[0.,v2[1]],[0.,v3[2]]]
atoms=asarray(zip(axs,ays,azs))

#Grid properties
nGridPoints=reduce(operator.mul,gridSize)

#Read in ELF data
elf=asarray(elf)
elf.shape=gridSize

#Neighbors
halfNeighbors=voronoiNeighbors(atoms=atoms,basis=basis,atypes=atypes,style='half')

print "Number of Neighbors before elimination:",sum([len(i) for i in halfNeighbors])

a=elf.shape
print "Average ELF value:",sum([sum([sum(line)/a[2] for line in plane])/a[1] for plane in elf])/a[0]

#Evaluate the ELF between each nieghbor pair
avgelfline,elflines,halfNeighbors=fieldNeighbors(atoms,atypes,basis,elf,gridSize,halfNeighbors)

#Cutoff neighbors that fall below the thresh-hold
Ninterp=len(avgelfline)
avgbig=zeros([len(bondCutoffs),Ninterp])
avgsmall=zeros([len(bondCutoffs),Ninterp])
nab=zeros(len(bondCutoffs))
nas=zeros(len(bondCutoffs))

cnt=0
for a,jNeighbors in enumerate(halfNeighbors):
    for b,jNeighb in enumerate(jNeighbors):
        vals=elflines[a][b]
        for j,bondCutoff in enumerate(bondCutoffs):
                bb=False
                if sum([1 for i in vals if i<bondCutoff])>0:
                    bb=True
                if bb:
                    nas[j]+=1
                    avgsmall[j]+=array(vals)                   
                else:
                    cnt+=1
                    nab[j]+=1
                    avgbig[j]+=array(vals)                   
for i in range(len(bondCutoffs)):
    if nas[i]==0:
        nas[i]+=1
    if nab[i]==0:
        nab[i]+=1
avgsmall=[avgsmall[i]/nas[i] for i in range(len(bondCutoffs))]
avgbig=[avgbig[i]/nab[i] for i in range(len(bondCutoffs))]
fige=pl.figure()
sml=fige.add_subplot(212)
big=fige.add_subplot(211)
print "Neighbors after elim:",cnt
for i,(avgs,avgb) in enumerate(zip(list(avgsmall),list(avgbig))):
    sml.plot(avgs,label=str(float(i)/10.))
    big.plot(avgb)

#pl.legend(["Good Neighbs","Bad Neighbs"])
sml.set_title("Average ELF value for non-Bonds")
sml.legend(title="cutoff",bbox_to_anchor=(1.10,1,))
big.set_title("Average ELF value for Bonds")

pl.show()
