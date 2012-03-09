#!/usr/bin/python

#This script rescales the desired POSCAR and then uploads its directory to gmice and starts a vasp run

import sys,os
#mine
import gmice_up,gmice_run

def usage():
    print "%s <local_dir> <remote_dir> <numprocs>"%sys.argv[0]
    exit(0)

if len(sys.argv)!=4:
    usage()

localdir = sys.argv[1].strip("/")+"/"
remotedir = sys.argv[2].strip("/")+"/"
numprocs = int(sys.argv[3])

#Copy to gmice
gmice_up.dir(localdir,remotedir)

#Run on gmice
gmice_run.vasp(remotedir+localdir,localdir[:15].strip("/"),numprocs)
