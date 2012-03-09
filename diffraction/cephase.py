#!/usr/bin/python
import sys
import pylab as pl
import os,glob
import pickle
from numpy import *
from math import asin
from scipy.signal import cspline1d,cspline1d_eval
from scipy.interpolate import UnivariateSpline
from scipy.optimize import leastsq
from optparse import OptionParser

#cephase.py -Adam Cadien, August 2010 GMU Computational Sciences and Informatics
#This program reads in a chi file, the settings for that chi file and
#a desired phase for cerium and then plots it over the chi file using
#the temperature, pressure and wavelength information for the chi file.


####################
#Arguement parser, should be the first thing called
def get_opts():
    usage="usage: %prog [options]"
    parser=OptionParser(usage=usage)
    parser.add_option("-i",dest="chifile",help="Chi file to plot in background.")
    parser.add_option("-p",dest="phase",help="Cerium phase to plot diffraction data of.")
    parser.add_option("-s",dest="settingsfile",help="Settings file, contains data for the chi file like the pressure, temperature, wavelength etc")
    
    (opts,args)=parser.parse_args()

    if opts.chifile is None:
        parser.print_help()
        print "Need a valid chi file."
        exit(0)

    if opts.phase is None:
        parser.print_help()
        print "Need a phase of cerium to plot."
        exit(0)
        
    if opts.settingsfile is None:
        parser.print_help()
        print "Need a valid settings file. Format is:"
        print "<chifilename>\\t<Temp>\\t<Lamda>\\t<Pressure>\\t<Ignored>"
        exit(0)

    return opts

###############
#Reads in a chi file, returns xxs and yys
def readchi(chifilename):
    f=open(chifilename,'r')
    i=0
    for line in f:
        if i<3:
            i+=1
            continue
        line=line.lstrip().rstrip()
        if i==3:
            i+=1
            length=int(line)
            break
    xxs=zeros(length)
    yys=zeros(length)
    i=0
    for line in f:
        line=line.lstrip().rstrip()
        [xxs[i],NULL,yys[i]]=line.split(' ')
        i+=1
    return (xxs,yys)

###############
#Reads in the settings file, returns temperature, pressure and wavelength
def readsettings(chifile,settingsfilename):
    f=open(settingsfilename,'r')
    found=False
    for line in f:
        if line[0]=='#':
            continue
        line=line.lstrip().rstrip()
        [fname,temperature,lamda,pressure,line]=line.split(None,4)
        if fname==chifile:
            return (float(temperature),float(pressure),float(lamda))
    print "Error: Unable to find "+chifile+" in settings file."
    exit(0)
          
##################
#Main Starts Here
##################
opts=get_opts()
chifile=opts.chifile
(xxs,yys)=readchi(chifile)
(T,P,L)=readsettings(chifile,opts.settingsfile)



