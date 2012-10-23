#!/usr/bin/python

#This script rescales the desired POSCAR and then uploads its directory to gmice and starts a vasp run

import sys,os
#mine
import rescaleposcar,gmice_up,gmice_run

def usage():
    print "%s <scale> <local_dir> <remote_dir> <numprocs>"%sys.argv[0]
    exit(0)

if len(sys.argv)!=5:
    usage()

perc = float(sys.argv[1])
localdir = sys.argv[2].strip("/")+"/"
remotedir = sys.argv[3].strip("/")+"/"
numprocs = int(sys.argv[4])

#Rescale
rescaleposcar.rescaleposcar(localdir+"POSCAR",perc)

#Copy 2 gmice
gmice_up.dir(localdir,remotedir)

#Run on gmice
gmice_run.vasp(remotedir+localdir,localdir[:15].strip("/"),numprocs)
