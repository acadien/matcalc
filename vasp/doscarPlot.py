#!/usr/bin/python

import sys
import subprocess

#for saving figures remotely
#import matplotlib
#matplotlib.use('Agg')

import pylab as pl
from scipy import *

def usage():
    print "doscarPlot.py <DOSCAR>"

if len(sys.argv)<2:
    usage()
    exit(0)

doscar=open(sys.argv[1]).readlines()

efermi=float(doscar[5].split()[-2])
Natoms=int(doscar.pop(0).split()[0])
doscar=doscar[4:]
NDOS=int(doscar.pop(0).split()[2])

#####################
# Get the total DOS #
#####################

#Check if spin was enabled for this calculation
spin=False
if len(doscar[0].split())==5:
    spin=True

tDOSenergy=list()
if spin:
    tDOSU=list()
    tDOSD=list()
    tDOSintegU=list()
    tDOSintegD=list()

    for i in range(NDOS):
        e,ud,dd,ui,di=map(float,doscar.pop(0).split())
        tDOSenergy.append(e)
        tDOSU.append(ud)
        tDOSD.append(dd)
        tDOSintegU.append(ui)        
        tDOSintegD.append(di)

else:
    tDOS=list()
    tDOSinteg=list()

    for i in range(NDOS):
        e,d,i=map(float,doscar.pop(0).split())
        tDOSenergy.append(e)
        tDOS.append(d)
        tDOSinteg.append(i)

###########################
# Get the per-orbital-DOS #
###########################

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

if len(poDOS) in [16,32]:#s,p,d,f
    orbs=[[0,1],[1,4],[4,9],[9,16]]
    labels=['s','p','d','f','total','integ']
elif len(poDOS) in [9,18]:#s,p,d
    orbs=[[0,1],[1,4],[4,9]]
    labels=['s','p','d','total','integ']
elif len(poDOS) in [4,8]:#s,p
    orbs=[[0,1],[1,4]]
    labels=['s','p','total','integ']
else:#s
    orbs=[[0,1]]
    labels=['s','total','integ']
colors=['blue','green','purple','red','black','gray']

#Shift everything over by the fermi Energy
poDOSenergy=[i-efermi for i in poDOSenergy]
tDOSenergy=[i-efermi for i in tDOSenergy]

pl.figure()
if spin:
    for i in range(len(orbs)):
        for io,o in enumerate(range(orbs[i][0],orbs[i][1])):
            if io==0:
                pl.plot(poDOSenergy,poDOS[o*2],c=colors[i],label=labels[i])
                pl.plot(poDOSenergy,poDOS[o*2+1]*-1,c=colors[i])
            else:
                pl.plot(poDOSenergy,poDOS[o*2],c=colors[i])
                pl.plot(poDOSenergy,poDOS[o*2+1]*-1,c=colors[i])

else:
    for i in range(len(orbs)):
        for io,o in enumerate(range(orbs[i][0],orbs[i][1])):
            if io==0:
                pl.plot(poDOSenergy,poDOS[o],c=colors[i],label=labels[i])
            else:
                pl.plot(poDOSenergy,poDOS[o],c=colors[i])

if spin:
    pl.plot(tDOSenergy,tDOSU,c=colors[-2],label=labels[-2])
    pl.plot(tDOSenergy,[x+y for x,y in zip(tDOSU,tDOSD)],c="orange",label="Sum of U&D",ls='--')
    #pl.plot(tDOSenergy,tDOSintegU,c=colors[-2],label=labels[-1],ls='-.')
    pl.plot(tDOSenergy,array(tDOSD)*-1,c=colors[-2])
    #pl.plot(tDOSenergy,array(tDOSintegD)*-1,c=colors[-2],ls='-.')
    
else:
    pl.plot(tDOSenergy,tDOS,c=colors[-2],label=labels[-2])
    #pl.plot(tDOSenergy,tDOSinteg,c=colors[-1],label=labels[-1],ls='-.')
pl.xlabel("Energy (eV)")
pl.ylabel("DOS")
pl.legend(loc=0)

mn=min(min(poDOSenergy),min(tDOSenergy))
mx=max(max(poDOSenergy),max(tDOSenergy))
pl.xlim([mn,mx])

pl.text(mn+(mx-mn)*0.01,max(tDOS)*0.98,"E-Fermi=%4.4feV"%efermi)
pl.plot([0,0],[0,max(tDOS)],c='black',ls=':',lw=2)

pl.title(sys.argv[1])

pl.show()



