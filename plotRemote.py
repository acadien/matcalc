#!/usr/bin/python

import os

#This script detects if a session is being run remotely:
#if remote - save figure to directory BASEDIR
#if local - plot it
#How to use this script: import *before* importing pylab or mpl
#use prshow(fname) instead of pylab.show()

REMOTESESSION_BASEDIR="~/Dropbox/"

try:
    os.environ['SSH_CLIENT']
except KeyError:
    #Local connection
    REMOTESESSION=False
else:
    REMOTESESSION=True
    import matplotlib
    matplotlib.use("Agg")


def prshow(fname="latestplot.png"):
    if REMOTESESSION:
        matplotlib.pyplot.save(REMOTESESSION_BASEDIR+fname)
    else:
        matplotlib.pyplot.show()
        
