#!/usr/bin/python

import sys
import pylab as pl
from scipy import *

def usage():
    print "doscarPlot.py <DOSCAR>"

if len(sys.argv)<2:
    usage()
    exit(0)

doscar=open(sys.argv[1]).readlines()

Natoms=int(doscar.pop(0).split()[0])
doscar=doscar[4:]
NDOS=int(doscar.pop(0).split()[2])

#Get the total DOS
tDOSenergy=list()
tDOS=list()
tDOSinteg=list()
for i in range(NDOS):
    e,d,i=map(float,doscar.pop(0).split())
    tDOSenergy.append(e)
    tDOS.append(d)
    tDOSinteg.append(i)

#Get the per-orbital-DOS
Norbs=len(doscar[2].split())-1
poDOSenergy=zeros(NDOS)
poDOS=zeros([NDOS,Norbs])
for a in range(Natoms):
    doscar.pop(0)
    for i in range(NDOS):
        line=map(float,doscar.pop(0).split())
        poDOSenergy[i]+=line[0]
        poDOS[i]+=line[1:]
    poDOS[i]/=Natoms
poDOSenergy/=Natoms
poDOS=poDOS.T

orbs=[[0,1],[1,4],[4,9],[9,16]]
labels=['s','p','d','f','total','integ']
colors=['blue','green','purple','red','black','black']
pl.figure()
for i in range(len(orbs)):
    for io,o in enumerate(range(orbs[i][0],orbs[i][1])):
        if io==0:
            pl.plot(poDOSenergy,poDOS[o],c=colors[i],label=labels[i])
        else:
            pl.plot(poDOSenergy,poDOS[o],c=colors[i])
pl.plot(tDOSenergy,tDOS,c=colors[-2],label=labels[-2])
pl.plot(tDOSenergy,tDOSinteg,c=colors[-1],label=labels[-1],ls='--')
pl.xlabel("Energy")
pl.ylabel("DOS")
pl.legend()
pl.xlim([min(min(poDOSenergy),min(tDOSenergy)),max(max(poDOSenergy),max(tDOSenergy))+8])
pl.title(sys.argv[1])
pl.show()



