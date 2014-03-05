#!/usr/bin/python

import sys
import subprocess

from scipy import *

def usage():
    print "%s <list of doscars>"%sys.argv[0].split("/")[-1]
    print "Writes files OUTCAR.col in respective directories"

if len(sys.argv)<2:
    usage()
    exit(0)

doscarFiles = sys.argv[1:]
outcarFiles = ["/".join(i.split("/")[:-1]+["OUTCAR"]) for i in doscarFiles]
TOoutputFiles = ["/".join(i.split("/")[:-1]+["TODOS.col"]) for i in doscarFiles]

for dcf,ocf,TOoutput in zip(doscarFiles,outcarFiles,TOoutputFiles):
    doscar=open(dcf,"r").readlines()

    #fermi energy
    efermi=float(subprocess.check_output("tail -n 2000 %s | grep E\-fermi"%ocf,shell=True).split("\n")[-2].split()[2])
    enableFermi=True

    #nAtoms, nDOS
    nAtom=int(doscar.pop(0).split()[0])
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
            tDOSU.append(ud/nAtom)
            tDOSD.append(dd/nAtom)
            tDOSintegU.append(ui/nAtom)        
            tDOSintegD.append(di/nAtom)

    else:
        tDOS=list()
        tDOSinteg=list()

        for i in range(NDOS):
            e,d,i=map(float,doscar.pop(0).split())
            tDOSenergy.append(e)
            tDOS.append(d/nAtom)
            tDOSinteg.append(i/nAtom)

    #Shift everything over by the fermi Energy
    tDOSenergy=[i-efermi for i in tDOSenergy]

    #TODOS.columns
    data=["Energy(eV) n(E)\n"]
    if spin:
        data+=["%lf %lf %lf\n"%(tDOSenergy[i],tDOSU[i],tDOSD[i]) for i in range(len(tDOSenergy))]
    else:
        data+=["%lf %lf\n"%(tDOSenergy[i],tDOS[i]) for i in range(len(tDOSenergy))]
    open(TOoutput,"w").writelines(data)

