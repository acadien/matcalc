#!/usr/bin/python

import sys,os

#This script copies the INCAR to INCAR_old and inserts a new INCAR for generating the DOSCAR from the WAVECAR (which must exist in the current directory).  It then submits this job to the queue to run on 1 processor.

def usage():
    print "repeatDOSCAR.py <dir with INCAR/WAVECAR etc>"

if len(sys.argv)!=2:
    usage()

wdir=sys.argv[1].strip("/")+"/"

os.system("cp %sINCAR %sINCAR_old"%(wdir,wdir))
os.system("cp /home/acadien/scripts/INCAR_dos %sINCAR"%wdir)
absdir=os.getcwd().strip("/")+"/"+wdir
absdir="/".join(absdir.split("/")[2:])
name=wdir[:15].strip("/")
os.chdir("/home/acadien")
os.system("./vsprun.sh %s %s 1"%(absdir,name))
