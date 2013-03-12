#!/usr/bin/python

import os

#This script detects if a session is being run remotely:
#if remote - save figure to directory BASEDIR
#if local - plot it
#use prmshow(fname) instead of mlab.show()

from mayavi import mlab                                                    

REMOTESESSION_BASEDIR="/home/acadien/Dropbox/"
try:
    os.environ['SSH_CLIENT']
except KeyError:
    #Local connection
    REMOTESESSION=False
else:
    REMOTESESSION=True
    mlab.options.offscreen = True


def prmshow(fname="latestplot.png"):
    if REMOTESESSION:
        mlab.savefig(REMOTESESSION_BASEDIR+fname)                 
    else:
        mlab.show()

        
