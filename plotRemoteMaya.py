#!/usr/bin/python

import os

#This script detects if a session is being run remotely:
#if remote - save figure to directory BASEDIR
#if local - plot it
#How to use this script: import *before* importing pylab or mpl
#use prshow(fname) instead of pylab.show()

REMOTESESSION_BASEDIR="/home/acadien/Dropbox/"
try:
    os.environ['SSH_CLIENT']
except KeyError:
    #Local connection
    REMOTESESSION=False
else:
    REMOTESESSION=True
    from mayavi import mlab                                                    
    mlab.options.offscreen = True


def prmshow(fname="latestplot.png"):
    if REMOTESESSION:
        mlab.savefig(REMOTESESSION_BASEDIR+fname)                 
    else:
        mlab.show()

        
