#!/usr/bin/python

import sys

#for saving figures remotely
#import matplotlib
#matplotlib.use('Agg')

import pylab as pl
from scipy import *
import numpy
#mine
from gaussFunctions import gauss1D

def usage():
    print "doscarPlot.py <DOSCAR> <OUTCAR (Efermi)>"

def gaussSmooth(xs,yrough,xi,yi,sig):
    #ys=[y*gaussNorm1D(x,mu,sig) for x,y in zip(xs,yrough)]

    Nx=len(xs)
    Nxi=len(xi)

    for j in range(Nx):
        if yrough[j]==0.0: continue
        low=max(0,j*5-50)
        high=min(Nxi,j*5+50)
        yi[low:high]+=gauss1D(array(xi[low:high]),yrough[j],xs[j],sig)
 
if len(sys.argv)<2:
    usage()
    exit(0)

doscar=open(sys.argv[1]).readlines()

if len(sys.argv)==3:
    try:
        efermi=float([line for line in open(sys.argv[2]).readlines() if "E-fermi" in line][-1].split()[2])
    except IndexError:
        print "ERROR: Can't find Fermi energy (E-fermi line in OUTCAR), exiting."
        exit(0)
    else:
        print "EFermi:",efermi

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

print poDOS.shape
if len(poDOS) in [16,32]:
    orbs={'s':[0,1],'p':[1,4],'d':[4,9],'f':[9,16]}
elif len(poDOS) in [9,18]:
    orbs={'s':[0,1],'p':[1,4],'d':[4,9]}
elif len(poDOS) in [4,8]:
    orbs={'s':[0,1],'p':[1,4]}
elif len(poDOS) in [1,2]:
    orbs={'s':[0,1]}
else:
    print "Length of DOS is weird (not one of 1,4,9,16). exiting."
    exit(0)
labels=['s','p','d','f','total','integ']
colors=['blue','green','purple','red','black','gray']
    
pl.figure()
"""
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
    pl.plot(tDOSenergy,tDOSintegU,c=colors[-2],label=labels[-1],ls='-.')
    pl.plot(tDOSenergy,array(tDOSD)*-1,c=colors[-2])
    pl.plot(tDOSenergy,array(tDOSintegD)*-1,c=colors[-2],ls='-.')
    
else:
    pl.plot(tDOSenergy,tDOS,c=colors[-2],label=labels[-2])
    pl.plot(tDOSenergy,tDOSinteg,c=colors[-1],label=labels[-1],ls='-.')
pl.xlabel("Energy")
pl.ylabel("DOS")
pl.legend(loc=0)

pl.plot([efermi,efermi],[0,max(tDOSinteg)],c='black',ls=':')
pl.xlim([min(min(poDOSenergy),min(tDOSenergy)),max(max(poDOSenergy),max(tDOSenergy))+8])
pl.title(sys.argv[1])
"""
#exit(0)
Ninterp=len(poDOSenergy)*5
smthEn=[(i/float(Ninterp)*(poDOSenergy[-1]-poDOSenergy[0]))+poDOSenergy[0] for i in range(Ninterp)]
smthDOSu,smthDOSd=zip(*[[zeros(Ninterp),zeros(Ninterp)] for i in range(len(orbs))])
#smthDDOSu,smthDDOSd=zeros(Ninterp),zeros(Ninterp)
#smthPDOSu,smthPDOSd=zeros(Ninterp),zeros(Ninterp)
#smthSDOSu,smthSDOSd=zeros(Ninterp),zeros(Ninterp)
smthTotDOSU=zeros(Ninterp)
smthTotDOSD=zeros(Ninterp)
sigma=0.25
#Loop over -f- orbitals and sum up contributions with gaussian smoothing
for i,orb in enumerate(orbs.keys()):
    for io,o in enumerate(range(*orbs[orb])):
        gaussSmooth(poDOSenergy,map(float,poDOS[2*o]),smthEn,smthDOSu[i],sigma)
        gaussSmooth(poDOSenergy,map(float,poDOS[2*o+1]),smthEn,smthDOSd[i],sigma)
gaussSmooth(tDOSenergy,tDOSU,smthEn,smthTotDOSU,sigma)
gaussSmooth(tDOSenergy,tDOSD,smthEn,smthTotDOSD,sigma)
smthEn=[i-efermi for i in smthEn]

for i in range(len(orbs)):
    pl.plot(smthEn,smthDOSu[i],c=colors[i],label=orbs.keys()[i])
    pl.fill_between(smthEn,0,smthDOSu[i],color=colors[i],label="f-U",alpha=0.7)
    pl.fill_between(smthEn,0,smthDOSd[i]*-1,color=colors[i],label="f-D",alpha=0.5)    


pl.plot(smthEn,smthTotDOSU,color="black",label="Total")
pl.plot(smthEn,smthTotDOSD*-1,color="black")
#pl.plot([0,0],[max(smthDOSd)*-1,max(smthDOSu)],c='black',ls=':')
#pl.xlim([-5,max(poDOSenergy)-efermi])
pl.xlabel("Energy")
pl.ylabel("DOS")
pl.legend(loc=0)
pl.show()



