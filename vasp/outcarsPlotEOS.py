#!/usr/bin/python

import plotRemote as pr #mine 

import sys,os,shutil,itertools
import lammps #requires lammps.so file exists and environment be setup
              #info here: http://lammps.sandia.gov/doc/Section_python.html
import pylab as pl
from numpy import linspace
from matplotlib import colors
import matplotlib.pyplot as plt

#mine
import poscar2lmpcfg,poscarGrow,poscarVolume

try:
    os.remove("POSCAR_eos")
except OSError:
    pass

bars2GPa=1./10000.
kB2GPa=1./10.
aa3GPa2eV=1./1602.17656

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
    return enrg,natoms,vol,pres

#Prepare LAMMPS single point energy calculations
def initLammpsCmds(potential):
    lmpelems=[]
    masses=[]

    if "adp" in potential:
        pairstyle="pair_style adp"
        lmpelems=open(potential,"r").readlines()[3].split()[1:]
    elif ("cbb" in potential) or ("CBB" in potential) or ("kawa" in potential):
        pairstyle="pair_style kawamura 7.0"
        masses=["mass 1 28.0855","mass 2 15.9994"]
    else:
        print "Error: Potential style not recognized"
        exit(0)

    preLammpsCommands=[\
        "units metal","boundary p p p","atom_style atomic", \
        "neighbor 0.3 bin","neigh_modify delay 5","box tilt large"]
    postLammpsCommands=masses+[\
        pairstyle, "pair_coeff * * %s %s"%(potential," ".join(lmpelems)), \
        "variable v equal vol",\
        #"minimize 1E-4 0.0 1000 10000",\
        "minimize 0 0 0 0",\
        ]
    return preLammpsCommands,postLammpsCommands

def lammpsGenerateE(vaspPOSCAR,preCmd,postCmd,vRatio):
    #Convert to lammps input configure format, rescale accordingly
    lammpsConfig="lmp.config"
    lammpsPOSCAR="POSCAR_eos"
    shutil.copyfile(vaspPOSCAR,lammpsPOSCAR)
    poscarGrow.poscarGrow(lammpsPOSCAR,lammpsPOSCAR,3,3,3)
    poscarVolume.ratio(lammpsPOSCAR,vRatio)
    poscar2lmpcfg.poscar2dump(lammpsPOSCAR,lammpsConfig)
    
    #Run lammps
    lmp=lammps.lammps()
    print "here"
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
    print sys.argv[0].split("/")[-1]+" <base eos directory> <optional:A B C scaling params> <optional: lammps potential>"
    print "Loops through all sub directories grabs OUTCARs and energies for EOS"
    print "If LAMMPS potential is given, compares LAMMPS and VASP EOS"

if len(sys.argv)<2:
    usage()
    exit(0)

basedir = sys.argv[1].rstrip("/")
lmppot=-1
A=1.0
B=1.0
C=0.0
if len(sys.argv)>2:
    if len(sys.argv)==3:
        lmppot=sys.argv[2]
    elif len(sys.argv)==5:
        A=float(sys.argv[2])
        B=float(sys.argv[3])
        C=float(sys.argv[4])
    elif len(sys.argv)==6:
        A=float(sys.argv[2])
        B=float(sys.argv[3])
        C=float(sys.argv[4])
        lmppot=sys.argv[5]

#Prepare LAMMPS single point energy calculations
phases=[i for i in os.listdir(basedir) if "eos" in i]

#VASP Data
ratios={}
Vvolumes={}
Vpressures={}
Venergies={}
Venthalpies={}
toremove=list()
for phase in phases:
    if len([i for i in os.listdir(basedir+"/"+phase)]) == 0:
        print phase,"directy empty, skipping and moving on"
        toremove.append(phase)
        continue

    ratios[phase]=[i for i in os.listdir(basedir+"/"+phase) if os.path.isdir(basedir+"/"+phase+"/"+i) and "original" not in i]

    es,natoms,vols,prss = zip(* \
        [outcarGrabFinalE("/".join([basedir,phase,rat,"OUTCAR"])) \
             for rat in ratios[phase]])
    natoms=natoms[0]
    Venergies[phase],Vvolumes[phase],Vpressures[phase] = \
        zip(* sorted(zip(es,vols,prss),key=lambda x:x[1]) )
    Venergies[phase] = [i*B-C for i in Venergies[phase]]
    Vvolumes[phase] = [i*A**3 for i in Vvolumes[phase]]
    Vpressures[phase] = [i*B/(A**3) for i in Vpressures[phase]]
    Venthalpies[phase] = [e + p*v*aa3GPa2eV for e,v,p in zip(Venergies[phase],Vvolumes[phase],Vpressures[phase])]

for i in toremove:
    phases.remove(i)

#LAMMPS Data
Lvolumes={}
Lpressures={}
Lenergies={}
Lenthalpies={}
if lmppot!=-1:
    preCmd,postCmd=initLammpsCmds(lmppot)
    for phase in phases:
        print phase
        #Volume ratios for lammps configurations
        minr=float(min(ratios[phase]))
        maxr=float(max(ratios[phase]))
        numr=30

        epv = [lammpsGenerateE("/".join([basedir,phase,"1.00/POSCAR"]),preCmd,postCmd,r) \
                   for r in [0.80,0.85,0.90,0.92,0.95,0.98,1.00,1.02,1.04,1.06,1.10,1.15,1.20]]
                   #for r in linspace(minr,maxr,numr)]



        Lenergies[phase],Lpressures[phase],Lvolumes[phase] = \
            zip(* sorted(epv,key=lambda x:x[2]) )
        
        Lenthalpies[phase] = [e + p*v*aa3GPa2eV for e,v,p in zip(Lenergies[phase],Lvolumes[phase],Lpressures[phase])]

#Plotting 
mcolors=[pl.cm.spectral(i) for i in linspace(0,0.9,len(phases))]
markers=['o','v','s','p','*','h','D']

subs=140
if lmppot!=-1:
    subs=220

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

def dictScatterInset(xdict,ydict,items,xlim,ylim,xticks,yticks):
    a = pl.axes([.185, .55, .3, .3], axisbg='gray')
    cs=itertools.cycle(mcolors)
    marks=itertools.cycle(markers)
    for i in items:
        c=cs.next()
        m=marks.next()
        pl.scatter(xdict[i],ydict[i],marker=m,c=c,label=i,s=40)
        pl.plot(xdict[i],ydict[i],c=c)
        pl.setp(a,xlim=xlim,ylim=ylim,xticks=xticks,yticks=yticks)
#Plot

pl.figure(figsize=(20,12))
#pl.subplot(subs+1)

dictScatter(Vvolumes,Venergies,phases)
if lmppot!=-1: dictPlot(Lvolumes,Lenergies,phases,"-",1.5)
pl.xlabel("Volume ($\AA^3 / atom$)",size=17)
pl.ylabel("Energy ($eV / atom$)",size=17)
pl.legend(loc=0,fontsize=12)


#pl.subplot(subs+2)
pl.figure(figsize=(20,12))
dictScatter(Vpressures,Venergies,phases)
if lmppot!=-1: dictPlot(Lpressures,Lenergies,phases,"-",1.5)
pl.xlabel("Pressure ($GPa$)",size=17)
pl.ylabel("Energy ($eV / atom$)",size=17)
pl.ylim([-4.55,-3])
pl.xlim([-15,50])
pl.legend(loc=0,fontsize=12)
#dictScatterInset(Vpressures,Venergies,phases,xlim=[-20,20],ylim=[-8,-7.5], xticks=[-20,0,20], yticks=[-8,-7.5])
#pr.prshow("EOS_PE.png")
pr.prshow("EOS_PVol.png")
exit(0)
#pl.subplot(subs+3)
pl.figure()
dictScatter(Vpressures,Venthalpies,phases)
#pl.xlim([0,500])
#pl.ylim([-8.5,-2])
if lmppot!=-1: dictPlot(Lpressures,Lenthalpies,phases,"-",1.5)
pl.xlabel("Pressure ($GPa$)",size=17)
pl.ylabel("Enthalpy ($eV / atom$)",size=17)
pl.legend(loc=0,fontsize=12)
#dictScatterInset(Vpressures,Venergies,phases,xlim=[0,100],ylim=[-8,-7], xticks=[0,20,40,60,80,100], yticks=[-8,-7])
#pr.prshow("EOS_PH.png")

#pl.subplot(subs+4)
pl.figure()
dictScatter(Vvolumes,Vpressures,phases)
if lmppot!=-1: dictPlot(Lvolumes,Lpressures,phases,"-",1.5)
pl.xlabel("Volume ($\AA^3 / atom$)",size=17)
pl.ylabel("Pressure ($GPa$)",size=17)
pl.legend(loc=0,fontsize=12)
#pr.prshow("EOS_PVol.png")


plt.tight_layout(pad=0.1,h_pad=0.1,w_pad=0.1)
pr.prshow("outcarsEOS.png")
