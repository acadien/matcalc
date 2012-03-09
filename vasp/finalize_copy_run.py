#!/usr/bin/python

#Creates a new directory taking for the POSCAR the last OUTCAR, copies it up to gmice and generates the ELFCAR, DOSCAR & AECCARS for a run on a 2x2x2 set of KPOINTS.

import sys,os
#mine
import outcar2poscar,gmice_up,gmice_run

def usage():
    print "%s <old_dir> <remote_dir>"%sys.argv[0]
    exit(0)

if len(sys.argv)!=3:
    usage()

local="/home/acadien/Documents/ascwork/scripts/vasp_scripts/"

stage1=sys.argv[1].strip("/")
stage2=stage1+"_fin/"
stage1+="/"

remotedir=sys.argv[2].strip("/")+"/"

#Prepare the new simulation locally
os.system("mkdir %s"%stage2)
os.system("cp %sPOTCAR %s"%(stage1,stage2))

#For DOS calculations.
#os.system("cp %sINCAR_dos %sINCAR"%(local,stage2))
#os.system("cp %sKPOINTS_dos %sKPOINTS"%(local,stage2))

#For precision structure refinement. (useful for Enthalpy vs Pressure)
#=======
incar=open("%s/INCAR"%stage1,"r").readlines()
incarprecise=open("%s/INCAR"%stage2,"w")
for line in incar:
    #need to change encut ?
    if "PREC" in line:
        line="PREC    = HIGH"
    incarprecise.write(line)


#=======


os.system("cp %s/KPOINTS %s"%(stage1,stage2))

outcar2poscar.outcar2poscar(stage1,"%sPOSCAR"%stage2,-1)

#Upload 
gmice_up.dir(stage2,remotedir)

#Run!
gmice_run.vasp(remotedir+stage2,stage2[:15].strip("/"),8)
