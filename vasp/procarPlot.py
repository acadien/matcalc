#!/usr/bin/python

import sys
import pylab as pl
#mine
from procarIO import readPROCAR

def usage():
    print "procarPlot.py <PROCAR File> <Value=avg or kp> <Xaxis=band or energy>"

if len(sys.argv)!=4:
    usage()
    exit(0)
procarfile=sys.argv[1]
yvals=sys.argv[2]
if yvals not in ["avg","kp"]:
    usage()
    print "Select a valid value to plot."
    exit(0)
xaxis=sys.argv[3]
if xaxis not in ["band","energy"]:
    usage()
    print "Select a valid xaxis."
    exit(0)

olabels,kpoints,ws,energy,occupancy,avgen,avgoc=readPROCAR(procarfile)

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
            if xaxis=='band':
                pl.plot(avgoc[i],label=temp,color=c)
            elif xaxis=='energy':
                pl.plot(avgen,avgoc[i],label=temp,color=c)
        else:
            if xaxis=='band':
                pl.plot(avgoc[i],color=c)
            elif xaxis=='energy':
                pl.plot(avgen,avgoc[i],color=c)
        porb=orb
    if xaxis=='band':
        pl.xlabel("Band#")
    elif xaxis=='energy':
        pl.xlabel("Energy (in eV)")
    pl.ylabel("DOS")
    pl.legend()
    if xaxis=='band':
        pl.xlim([0,len(avgen)+200])
        pl.savefig("DOS_vsBand_from%s"%procarfile)
    elif xaxis=='energy':
        pl.xlim([min(avgen),max(avgen)+12.0])
        pl.savefig("DOS_vsEnergy_from%s"%procarfile)

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
                if xaxis=='band':
                    pl.plot(occupancy[kp][i],label=temp,color=c)
                elif xaxis=='energy':
                    pl.plot(energy[kp],occupancy[kp][i],label=temp,color=c)
            else:
                if xaxis=='band':
                    pl.plot(occupancy[kp][i],color=c)
                elif xaxis=='energy':
                    pl.plot(energy[kp],occupancy[kp][i],color=c)
            porb=orb
        if xaxis=='band':
            pl.xlabel("Band#")
        elif xaxis=='energy':
            pl.xlabel("Energy (in eV)")
        pl.ylabel("DOS")
        pl.title("KPoint %d at (%5.5g,%5.5g,%5.5g)"%(kp,kpoints[kp][0],kpoints[kp][1],kpoints[kp][2]));
        pl.legend()
        if xaxis=='band':
            pl.xlim([0,len(avgen)+200])
            pl.savefig("DOS_vsBand_kp%d_from%s"%(kp,procarfile))
        elif xaxis=='energy':
            pl.xlim([min(avgen),max(avgen)+12.0])
            pl.savefig("DOS_vsEnergy_kp%d_from%s"%(kp,procarfile))

