#!/usr/bin/python

import plotRemote as pr

import sys
import pylab as pl
#from numpy import *
#mine
import procarIO


def usage():
    print "procarPlot.py <PROCAR File> <Value=avg or kp> <Interp = gauss or point>"

if len(sys.argv)!=4:
    usage()
    exit(0)
procarfile=sys.argv[1]
yvals=sys.argv[2]
if yvals not in ["avg","kp"]:
    usage()
    print "Select a valid value to plot."
    exit(0)
interp=sys.argv[3]
if interp not in ["gauss","point"]:
    usage()
    print "Select a valid interpolation/summation method."
    exit(0)

if interp=='gauss':
    olabels,kpoints,ws,energy,occupancy,enGrid,ocGrid=procarIO.read(procarfile)
else:
    olabels,kpoints,ws,energy,occupancy,enGrid,ocGrid=procarIO.read(procarfile,sigma=0)

Norbs=len(olabels)
porb='q'
if yvals=='avg':

    for i in range(Norbs):
        orb=olabels[i][0]
        if orb=='s':
            c='blue'
        elif orb=='p':
            c='green'
        elif orb=='d':
            c='purple'
        elif orb=='f':
            c='red'
        elif orb=='t':
            c='black'
        if orb!=porb:
            if orb=='t': temp='total'
            else: temp=orb
            pl.plot(enGrid,ocGrid[i],label=temp,color=c)
        else:
            pl.plot(enGrid,ocGrid[i],color=c)
        porb=orb

    pl.xlabel("Energy (in eV)")
    pl.ylabel("DOS")
    pl.legend()
    pr.prshow("DOS_from%s"%procarfile)

elif yvals=='kp':
    for kp in range(len(kpoints)):
        pl.figure()
        for i in range(Norbs):
            orb=olabels[i][0]
            if orb=='s':
                c='blue'
            elif orb=='p':
                c='green'
            elif orb=='d':
                c='purple'
            elif orb=='f':
                c='red'
            elif orb=='t':
                c='black'
            if orb!=porb:
                if orb=='t': temp='total'
                else: temp=orb
                pl.plot(energy[kp],occupancy[kp,:,i],label=temp,color=c)
            else:
                pl.plot(energy[kp],occupancy[kp,:,i],color=c)
            porb=orb
        pl.xlabel("Energy (in eV)")
        pl.ylabel("DOS")
        pl.title("KPoint %d at (%5.5g,%5.5g,%5.5g)"%(kp,kpoints[kp][0],kpoints[kp][1],kpoints[kp][2]));
        pl.legend(loc=0)
            #pl.xlim([min(avgen),max(avgen)+12.0])
        pl.savefig("DOS_kp%d_from%s"%(kp,procarfile))

