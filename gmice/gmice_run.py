#!/usr/bin/python

import os

def vasp(remotedir,name,numprocs):
    os.system("ssh gmice.gmu.edu \"vsprun.sh %s %s %d\""%(remotedir,name,numprocs))
