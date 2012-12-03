#!/usr/bin/python
import sys,os

#mine
import plotRemote as pr
import pylab as pl
import itertools
#mine
import colors

def usage():
    print "Usage:"
    print sys.argv[0].split("/")[-1]+" <base eos directory>"
    print "Loops through all sub directories grabs OUTCARs and energies for EOS"

if len(sys.argv)<2:
    usage()
    exit(0)

def outcarGrabFinalE(outcar):
    natoms=0
    vol=0
    pres=0
    for line in open(outcar):
        if "NIONS" in line:
            natoms=int(line.split()[-1])
        if "volume of cell" in line:
            vol=float(line.split()[-1])
        if "external pressure" in line:
            pres=float(line.split()[3])/10.0 #convert kB to GPa
        if "free  energy" in line:
            return float(line.split()[-2]),natoms,vol,pres
    print "Error: simulation didn't finish (not final TOTEN) in file %s"%outcar
    exit(0)

basedir = sys.argv[1].rstrip("/")
phases=[i for i in os.listdir(basedir) if "eos" in i[-3:]]

ratios={}
volumes={}
pressures={}
energies={}
for phase in phases:
    ratios[phase]=[i for i in os.listdir(basedir+"/"+phase) if "." in i]
    es,natoms,vols,prss = zip(* \
        [outcarGrabFinalE(basedir+"/"+phase+"/"+rat+"/OUTCAR") \
             for rat in ratios[phase]])
    natoms=natoms[0]
    es,vols,prss=zip(*sorted(zip(es,vols,prss),key=lambda x:x[1]))
    energies[phase]=[i/natoms for i in es]
    volumes[phase]=[i/natoms for i in vols]
    pressures[phase]=prss #don't divide by # atoms

colors=[colors.float2rgb(i,0,len(phases)) for i in range(len(phases))]
markers=['o','v','s','p','*','h','D']

#Plotting E vs Vol
pl.figure()
cs=itertools.cycle(colors)
marks=itertools.cycle(markers)
for phase in phases:
    c=cs.next()
    m=marks.next()
    pl.plot(volumes[phase],energies[phase],label=phase,c=c,mfc=c,marker=m)
pl.xlabel("Volume ($\AA / atom$)")
pl.legend(loc=0)
pl.ylabel("Energy ($eV / atom$)")
pr.prshow("EVolEOS.png")

#Plotting P vs Vol
pl.figure()
cs=itertools.cycle(colors)
marks=itertools.cycle(markers)
for phase in phases:
    c=cs.next()
    m=marks.next()
    pl.plot(volumes[phase],pressures[phase],label=phase,c=c,mfc=c,marker=m)
pl.xlabel("Volume ($\AA / atom$)")
pl.legend(loc=0)
pl.ylabel("Pressure ($GPa$)")
pr.prshow("PVolEOS.png")
 
#pl.savefig("/home/acadien/Dropbox/EVolEOS.png")
#pl.show()
