#!/usr/bin/python

import sys,os,re

def usage():
    print "%s <heating_dir> <cooling_dir>"%sys.argv[0].split("/")[-1]
    print "Copies KPOINTS/CONTCAR/POTCAR/INCAR from heating_dir to cooling_dir, reverses temperatures in INCAR"

if len(sys.argv)!=3:
    usage()
    exit(0)

heatdir=sys.argv[1].rstrip("/")+"/"
cooldir=sys.argv[2].rstrip("/")+"/"

if os.path.exists(cooldir):
    print "Please delete or rename cooling directory before running."
    exit(0)
if not os.path.exists(heatdir):
    print "Unable to locate heating directory."
    exit(0)

os.mkdir(cooldir)

comm=['cp',heatdir+"INCAR",cooldir+"INCAR"]
os.spawnvpe(os.P_WAIT, 'cp', comm, os.environ)
comm=['cp',heatdir+"CONTCAR",cooldir+"POSCAR"]
os.spawnvpe(os.P_WAIT, 'cp', comm, os.environ)
comm=['cp',heatdir+"KPOINTS",cooldir+"KPOINTS"]
os.spawnvpe(os.P_WAIT, 'cp', comm, os.environ)
comm=['cp',heatdir+"POTCAR",cooldir+"POTCAR"]
os.spawnvpe(os.P_WAIT, 'cp', comm, os.environ)

#Grab beginning and ending temperatures and swap them
incar=open(cooldir+"INCAR","r").readlines()
begT=-1
endT=-1
begLine=-1
endLine=-1
for i in range(len(incar)):
    if "TEBEG" in incar[i] and "TEEND" in incar[i]:
        begLine=endLine=i
        begT=re.findall(".*TEBEG\s*=\s*(\d+\.*\d+).*",incar[i])[0]
        endT=re.findall(".*TEEND\s*=\s*(\d+\.*\d+).*",incar[i])[0]
        break
    if "TEBEG" in incar[i]:
        begLine=i
        begT=re.findall(".*TEBEG\s*=\s*(\d+\.*\d+).*",incar[i])[0]
    if "TEEND" in incar[i]:
        endLine=i
        endT=re.findall(".*TEEND\s*=\s*(\d+\.*\d+).*",incar[i])[0]
    if begLine!=-1 and endLine!=-1:
        break

if begLine==endLine:
    incar[begLine]="TEBEG = %s ; TEEND = %s\n"%(endT,begT)
else:
    incar[begLine]="TEBEG = %s\n"%endT
    incar[endLine]="TEEND = %s\n"%begT
open(cooldir+"INCAR","w").writelines(incar)
