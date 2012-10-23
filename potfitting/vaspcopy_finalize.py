#!/usr/bin/python

import sys,os,re

def usage():
    print "%s <init_dir> <final_dir> <final-kpoint-file>"%sys.argv[0].split("/")[-1]
    print "Copies POSCAR/POTCAR/INCAR from init_dir to final_dir, sets precision to HIGH in INCAR, and inserts new kpoint file."

if len(sys.argv)!=4:
    usage()
    exit(0)

initdir=sys.argv[1].rstrip("/")+"/"
finaldir=sys.argv[2].rstrip("/")+"/"
kpointfil=sys.argv[3]

if os.path.exists(finaldir):
    print "Please delete or rename final directory before running."
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
comm=['cp',initdir+"POSCAR",finaldir+"POSCAR"]
os.spawnvpe(os.P_WAIT, 'cp', comm, os.environ)
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
