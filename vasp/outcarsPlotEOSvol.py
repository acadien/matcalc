#!/usr/bin/python
import sys,os

#for saving figures remotely
import matplotlib
matplotlib.use('Agg')

import pylab as pl

if len(sys.argv)<2:
    print "Usage:"
    print sys.argv[0].split("/")[-1]+" <base eos directory>"
    print "Loops through all sub directories grabs OUTCARs and energies for EOS"

def outcarGrabFinalE(outcar):
    natoms=0
    vol=0
    for line in open(outcar):
        if "NIONS" in line:
            natoms=int(line.split()[-1])
        if "volume of cell" in line:
            vol=float(line.split()[-1])
        if "free  energy" in line:
            return float(line.split()[-2]),natoms,vol
    print "Error: simulation didn't finish (not final TOTEN) in file %s"%outcar
    exit(0)

basedir = sys.argv[1].rstrip("/")
phases=[i for i in os.listdir(basedir) if "eos" in i[-3:]]

ratios={}
volumes={}
energies={}
for phase in phases:
    ratios[phase]=[i for i in os.listdir(basedir+"/"+phase) if "." in i]
    es,natoms,vols = zip(* \
        [outcarGrabFinalE(basedir+"/"+phase+"/"+rat+"/OUTCAR") \
             for rat in ratios[phase]])
    natoms=natoms[0]
    es,vols=zip(*sorted(zip(es,vols),key=lambda x:x[1]))
    energies[phase]=[i/natoms for i in es]
    volumes[phase]=[i/natoms for i in vols]

pl.figure()
for phase in phases:
    pl.plot(volumes[phase],energies[phase],label=phase)
pl.xlabel("Volume ($\AA / atoms$)")
pl.ylabel("Energy ($eV / atoms$)")
pl.legend(loc=0)
pl.savefig("/home/acadien/Dropbox/EVolEOS.png")
#pl.show()
