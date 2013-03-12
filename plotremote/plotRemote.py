#!/usr/bin/python

import os
import matplotlib

#This script detects if a session is being run remotely:
#if remote - save figure to directory BASEDIR
#if local - plot it
#How to use this script: import *before* importing pylab or mpl
#use prshow(fname) instead of pylab.show()

REMOTESESSION_BASEDIR = "/home/acadien/Dropbox/"

try:
    os.environ['SSH_CLIENT']
except KeyError:
    #Local connection
    REMOTESESSION = False
else:
    REMOTESESSION = True
    import matplotlib
    matplotlib.use("Agg")

import socket
if socket.gethostname()=="mozart":
    REMOTESESSION = False
    import matplotlib
    matplotlib.use('WX')

def prshow(fname="latestplot.png"):
    if REMOTESESSION:
        matplotlib.pyplot.savefig(REMOTESESSION_BASEDIR + fname)
    else:
        matplotlib.pyplot.show()
