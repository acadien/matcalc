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
#mine
from orderParam import coordinationNumber,bondOrientation
import poscarIO

orderParams={"CN":coordinationNumber, \
             "BO":bondOrientation    }

def usage():
    print "%s <POSCAR file> <order parameter>"%sys.argv[0].split("/")[-1]
    print "Order Parameter must be one of:"
    print "   CN : Coordination Number"
    print "   BO# : Bond Orientation (Q) with l=#"
    print ""

if len(sys.argv) < 3:
    usage()
    exit(0)

opFlag = False
lval=0
if len(sys.argv)==3:
    op = sys.argv[2]
    if op[:2]=="BO":
        lval=int(op[-1])
        op="BO"
    opFlag = True

poscar=open(sys.argv[1],"r").readlines()

[basis,atypes,atoms,head,poscar] = poscarIO.read(poscar)
ax,ay,az=zip(*atoms)
v1,v2,v3=basis
j=0
types=list()
for i in atypes:
    types+=[j+1]*i
    j+=1

#Get the order parameter
ops = orderParams[op](array(atoms),array(basis),lval)

if op=="BO":
    pl.hist(ops,bins=int(sqrt(len(ops))),normed=True,histtype='step')
    pl.xlabel(r"Bond Order $Q_%d$"%lval)
    pl.ylabel(r"$P ( Q_%d )$"%lval)
elif op=="CN":
    pl.hist(ops,bins=range(min(ops),max(ops)+2),normed=True,histtype='bar',align='left',rwidth=0.8)
    pl.xticks(range(min(ops),max(ops)+1))
    pl.xlabel(r"Coordination Number")
    pl.ylabel(r"P(CN)")

pr.prshow("%sOrderParam.png"%sys.argv[1])
