#!/usr/bin/python

import sys,os,re
#mine
import outcar2poscar

def usage():
    print "%s <init_dir> <final_dir> <outcar configuration> <final-kpoint-file>"%sys.argv[0].split("/")[-1]
    print "Copies POTCAR/INCAR from init_dir to final_dir, sets precision to HIGH in INCAR, inserts new kpoint file, and grabs a new POSCAR from the OUTCAR."

if len(sys.argv)!=5:
    usage()
    exit(0)

initdir=sys.argv[1].rstrip("/")+"/"
finaldir=sys.argv[2].rstrip("/")+"/"
configNum=int(sys.argv[3])
kpointfil=sys.argv[4]

if os.path.exists(finaldir):
    print "Delete or rename final directory before running."
    exit(0)
if not os.path.exists(initdir):
    print "Error: Unable to locate initial directory."
    exit(0)
try:
    open(kpointfil) 
except IOError:
    print "Error: Unable to locate initial KPOINT file."
    exit(0)

os.mkdir(finaldir)

comm=['cp',initdir+"INCAR",finaldir+"INCAR"]
os.spawnvpe(os.P_WAIT, 'cp', comm, os.environ)
#comm=['cp',initdir+"POSCAR",finaldir+"POSCAR"]
#os.spawnvpe(os.P_WAIT, 'cp', comm, os.environ)
comm=['cp',initdir+"POTCAR",finaldir+"POTCAR"]
os.spawnvpe(os.P_WAIT, 'cp', comm, os.environ)

comm=['cp',kpointfil,finaldir+"KPOINTS"]
os.spawnvpe(os.P_WAIT, 'cp', comm, os.environ)

#Change INCAR to highest accuracy and set NSW to 0
incar=open(finaldir+"INCAR","r").readlines()
for i in range(len(incar)):
    if "PREC" in incar[i]:
        incar[i]="PREC = high\n"
    if "IBRION" in incar[i]:
        incar[i]="IBRION = -1\n"
    if "NSW" in incar[i]:
        incar[i]="NSW = 0\n"
open(finaldir+"INCAR","w").writelines(incar)

#Make the new POSCAR from the selected configuration in the OUTCAR
outcar2poscar.outcar2poscar(initdir+"OUTCAR",finaldir+"POSCAR",configNum)
