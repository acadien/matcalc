#!/usr/bin/python

import plotRemote as pr #mine 

import sys,os,shutil,itertools
import lammps #requires lammps.so file exists and environment be setup
              #info here: http://lammps.sandia.gov/doc/Section_python.html
import pylab as pl
from numpy import linspace
from matplotlib import colors

#mine
import poscar2lmpcfg,poscarGrow,poscarVolume

bars2GPa=1./10000.
kB2GPa=1./10.

#Grab thermodynamic values from VASP simulations
def outcarGrabFinalE(outcar):
    natoms=0
    vol=0
    pres=0
    enrg=0
    for line in open(outcar):
        if "NIONS" in line:
            natoms=int(line.split()[-1])
        if "volume of cell" in line:
            vol=float(line.split()[-1])/natoms
        if "free  energy" in line:
            enrg=float(line.split()[-2])/natoms
        if "external pressure" in line:
            pres=float(line.split()[3])*kB2GPa
            return enrg,natoms,vol,pres
    print "Error: simulation didn't finish (not final TOTEN) in file %s"%outcar
    exit(0)

#Prepare LAMMPS single point energy calculations
def initLammpsCmds(potential):
    lmpelems=open(potential,"r").readlines()[3].split()[1:]
    preLammpsCommands=[\
        "units metal","boundary p p p","atom_style atomic", \
        "neighbor 0.3 bin","neigh_modify delay 5"]
    postLammpsCommands=[\
        "pair_style adp", "pair_coeff * * %s %s"%(potential," ".join(lmpelems)), \
        "variable v equal vol",\
        "minimize 0.0 0.0 0 0",\
        ]
    return preLammpsCommands,postLammpsCommands

def lammpsGenerateE(vaspPOSCAR,preCmd,postCmd,vRatio):
    #Convert to lammps input configure format, rescale accordingly
    lammpsConfig="lmp.config"
    lammpsPOSCAR="POSCAR_eos"
    shutil.copyfile(vaspPOSCAR,lammpsPOSCAR)
    poscarGrow.poscarGrow(lammpsPOSCAR,lammpsPOSCAR,2,1,1)
    poscarVolume.ratio(lammpsPOSCAR,vRatio)
    poscar2lmpcfg.poscar2dump(lammpsPOSCAR,lammpsConfig)
    
    #Run lammps
    lmp=lammps.lammps()
    map(lambda x:lmp.command(x),preCmd)
    lmp.command("read_data %s"%lammpsConfig)
    map(lambda x:lmp.command(x),postCmd)
    
    #Grab data from LAMMPS
    natom = lmp.get_natoms()
    pe = lmp.extract_compute("thermo_pe",0,0)/natom 
    prs = lmp.extract_compute("thermo_press",0,0)*bars2GPa 
    vol = lmp.extract_variable("v",0,0)/natom
    
    os.remove(lammpsConfig)
    os.remove(lammpsPOSCAR)
    
    return pe,prs,vol

#Processing input
def usage():
    print "Usage:"
    print sys.argv[0].split("/")[-1]+" <base eos directory> <optional: lammps potential>"
    print "Loops through all sub directories grabs OUTCARs and energies for EOS"
    print "If LAMMPS potential is given, compares LAMMPS and VASP EOS"

if len(sys.argv)<2:
    usage()
    exit(0)

basedir = sys.argv[1].rstrip("/")
lmppot=-1
if len(sys.argv)>2:
    lmppot=sys.argv[2]

#Prepare LAMMPS single point energy calculations
phases=[i for i in os.listdir(basedir) if "eos" in i[-3:]]

#VASP Data
ratios={}
Vvolumes={}
Vpressures={}
Venergies={}
for phase in phases:
    ratios[phase]=[i for i in os.listdir(basedir+"/"+phase) if "." in i]
    es,natoms,vols,prss = zip(* \
        [outcarGrabFinalE("/".join([basedir,phase,rat,"OUTCAR"])) \
             for rat in ratios[phase]])
    natoms=natoms[0]
    Venergies[phase],Vvolumes[phase],Vpressures[phase] = \
        zip(* sorted(zip(es,vols,prss),key=lambda x:x[1]) )

#LAMMPS Data
Lvolumes={}
Lpressures={}
Lenergies={}
if lmppot!=-1:
    preCmd,postCmd=initLammpsCmds(lmppot)
    for phase in phases:

        #Volume ratios for lammps configurations
        minr=float(min(ratios[phase]))
        maxr=float(max(ratios[phase]))
        numr=30

        epv = [lammpsGenerateE("/".join([basedir,phase,"1.00/POSCAR"]),preCmd,postCmd,r) \
                   for r in linspace(minr,maxr,numr)]

        Lenergies[phase],Lpressures[phase],Lvolumes[phase] = \
            zip(* sorted(epv,key=lambda x:x[2]) )

#Plotting 
mcolors=[pl.cm.spectral(i) for i in linspace(0,0.9,len(phases))]
markers=['o','v','s','p','*','h','D']

subs=130
#if lmppot!=-1:
#    subs=220

def dictPlot(xdict,ydict,items,ls,lw):
    cs=itertools.cycle(mcolors)
    marks=itertools.cycle(markers)
    for i in items:
        c=cs.next()
        m=marks.next()
        pl.plot(xdict[i],ydict[i],c=c,ls=ls,lw=lw)#mfc=c,marker=m,ls=ls,lw=lw)


def dictScatter(xdict,ydict,items):
    cs=itertools.cycle(mcolors)
    marks=itertools.cycle(markers)
    for i in items:
        c=cs.next()
        m=marks.next()
        pl.scatter(xdict[i],ydict[i],marker=m,c=c,label=i,s=40)
        pl.plot(xdict[i],ydict[i],marker=m,c=c,lw=1.5)

#VASP Plot
pl.subplot(subs+1)
dictScatter(Vvolumes,Venergies,phases)
if lmppot!=-1: dictPlot(Lvolumes,Lenergies,phases,"-",1.5)
pl.xlabel("Volume ($\AA / atom$)")
pl.ylabel("Energy ($eV / atom$)")
pl.legend(loc=0,fontsize=10)

pl.subplot(subs+2)
dictScatter(Vvolumes,Vpressures,phases)
if lmppot!=-1: dictPlot(Lvolumes,Lpressures,phases,"-",1.5)
pl.xlabel("Volume ($\AA / atom$)")
pl.ylabel("Pressure ($GPa$)")

pl.subplot(subs+3)
dictScatter(Vpressures,Venergies,phases)
if lmppot!=-1: dictPlot(Lpressures,Lenergies,phases,"-",1.5)
pl.xlabel("Pressure ($GPa$)")
pl.ylabel("Energy ($eV / atom$)")

"""
#VASP Plot
pl.subplot(subs+1)
dictPlot(Vvolumes,Venergies,phases)
pl.xlabel("Volume ($\AA / atom$)")
pl.ylabel("Energy ($eV / atom$)")
pl.legend(loc=0,fontsize=10)

pl.subplot(subs+2)
dictPlot(Vvolumes,Vpressures,phases)
pl.xlabel("Volume ($\AA / atom$)")
pl.ylabel("Pressure ($GPa$)")

#LAMMPS Plot
if lmppot!=-1:
    pl.subplot(subs+3)
    dictPlot(Lvolumes,Lenergies,phases)
    pl.xlabel("Volume ($\AA / atom$)")
    pl.ylabel("Energy ($eV / atom$)")

    pl.subplot(subs+4)
    dictPlot(Lvolumes,Lpressures,phases)
    pl.xlabel("Volume ($\AA / atom$)")
    pl.ylabel("Pressure ($GPa$)")

    pl.subplot(subs+5)
    dictPlotDiff(Vvolumes,Lvolumes,Venergies,Lenergies,phases)

    pl.subplot(subs+6)
    dictPlotDiff(Vvolumes,Lvolumes,Vpressures,Lpressures,phases)
"""
#pl.savefig("/home/acadien/Dropbox/EVolEOS.png")
#pl.show()
pr.prshow("PVolEOS.png")
