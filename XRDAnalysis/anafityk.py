#!/usr/bin/python
import pylab as pl
import os,glob,sys
from numpy import *
from optparse import OptionParser

##################
#func class
##################
class Func:
    def __init__(self):
        self.type=""
        self.center=-1.0
        self.height=-1.0
        self.area=-1.0
        self.FWHM=-1.0
        self.params=list()
    
    def set_props(self,vals):
        self.type=str(vals[1])
        self.center=double(vals[2])
        self.height=double(vals[3])
        self.area=double(vals[4])
        self.FWHM=double(vals[5])
        self.params=vals[6:]
            
    def set_props_quad(self,vals):
        self.type=str(vals[1])
        self.center=-1.0
        self.height=-1.0
        self.area=-1.0
        self.FWHM=-1.0
        self.params=vals[5:]

##################
#Main Starts Here
##################
usage="usage: %prog [options]"
parser=OptionParser(usage=usage)
parser.add_option("-i",dest="filename",help="Peaks data file name, should be output from fityk program")

(opts,args)=parser.parse_args()

if opts.filename is None:
    parser.print_help()
    print "\nError: Need an input peaks file from -i"
    exit(0)
else:
    fykfilename=opts.filename

curves=list()

fykf=open(fykfilename,"r")
for line in fykf:
    if line[0]=='%':
        curves.append(Func())
        params = [i for i in line.rstrip().split(' ') if i!='']
        if params[1]=='Quadratic':
            curves[-1].set_props_quad(params)
        else:
            curves[-1].set_props(params)

###At this point the curves list contains all the curves in the peaks list, specifically all the peak centers. can compare these centers with MgO centers
