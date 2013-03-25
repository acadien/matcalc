#!/usr/bin/python

import sys
#mine
import plotRemote as pr

from colors import float2rgb
from scipy import array
from mayavi import mlab
from numpy import linspace
import matplotlib.cm as cm
import pylab as pl
from math import sqrt
import re
#mine
from orderParam import coordinationNumber,bondOrientation,radialDistribution,bondAngleCorr
import poscarIO

orderParams={"CN":  coordinationNumber, \
             "BO":  bondOrientation , \
             "RDF": radialDistribution , \
             "BA":  bondAngleCorr }

def usage():
    print "%s <order parameter> <POSCAR Files (space delim)>"%sys.argv[0].split("/")[-1]
    print "Order Parameter must be one of:"
    print "   CN  : Coordination Number"
    print "   BO# : Bond Orientation (Q) with l=#"
    print "   RDF : Radial Distribution Function"
    print "   BA# : Bond Angle Correlation"
    print ""

if len(sys.argv) < 3:
    usage()
    exit(0)

lval=0

op = sys.argv[1]
if op[:2] in ["BO","BA"]:
    lval=int(op[-1])
    op=op[:2]

orderVals=list()
poscarNames=sys.argv[2:]
#Sort POSCAR names only based on the numbers in them
try:
    poscarNumbers=map(lambda x:float("".join(re.findall('\d+',x))),poscarNames)
    if len(poscarNumbers) == len(poscarNames):
        poscarNames=zip(*sorted(zip(poscarNames,poscarNumbers),key=lambda x:x[1]))[0]
except ValueError:
    pass


#Gather the order parameters
for pn in poscarNames:
    poscar=open(pn,"r").readlines()
    [basis,atypes,atoms,head,poscar] = poscarIO.read(poscar)

    #Get the order parameter
    orderVals.append(orderParams[op](array(atoms),array(basis),lval,debug=False))

#======================================================
#                       Plot!
#======================================================
if op=="BO":
    for i,ov in enumerate(orderVals):
        vals,bins,dummy = pl.hist(ov,bins=int(sqrt(len(ov)))*5,normed=True,visible=False)#,histtype='step')
        pl.plot(array(bins[:-1]),array(vals)+0.25*i)
    pl.xlabel(r"Bond Order $Q_%d$"%lval)
    pl.ylabel(r"$P ( Q_%d )$"%lval)

elif op=="CN":
    pl.hist(orderVals,bins=range(0,16),normed=True,histtype='bar',align='left',rwidth=0.8)
    pl.legend(poscarNames,loc=0)
    pl.xticks(range(min(map(min,orderVals)),max(map(max,orderVals))+1))
    pl.xlabel(r"Coordination Number")
    pl.ylabel(r"P(CN)")

elif op=="RDF":
    for ov in orderVals:
        pl.plot(ov[0],ov[1])
    pl.xlabel(r"R $( \AA )$")
    pl.ylabel("# Bonds")

elif op=="BA":
    for ov in orderVals:
        pl.plot(ov[0],ov[1])
    pl.xlabel(r"R $( \AA )$")
    pl.ylabel(r"$G_%d ( r )$"%lval)

pr.prshow("%s_chart_%s.png"%(sys.argv[1],op))
