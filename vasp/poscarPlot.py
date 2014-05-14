#!/usr/bin/python

import sys
#mine
import plotRemoteMaya as prm

from mayavi import mlab
import numpy as np
import pylab as pl
from math import *
#mine
from orderParam import coordinationNumber,bondOrientation,tetrahedral
from outcarPlotMSDAtom import outcarMeanSquareDisplaceAtom
from datatools import windowAvg
import poscarIO

orderParams={"CN":coordinationNumber, \
             "BO":bondOrientation, \
             "TET":tetrahedral,}
#MSD is not an order parameter but useful for labeling atoms at high temperatures.
#FILE is not an order parameter but applies colors based data from a file

def usage():
    print "----------------------------------------------------------------------------"
    print "%s <order parameter> <POSCAR file> <OPT: OUTCAR, REF>"%sys.argv[0]
    print "----------------------------------------------------------------------------"
    print "Order Parameter can be one of:"
    print "     CN : Coordination Number"
    print "    BO# : Bond Orientation (Q) with l=#"
    print "    TET : Tetrahedral order parameter (Sg)"
    print "    MSD : Mean Square Displacement, requires OUTCAR (annealed MD simulation)"
    print "   FILE : Parses a file (column, row, or CSV data) and applies coloring from that data"
    print "----------------------------------------------------------------------------"
    print "Flags can be used anywhere in args:"
    print "-rectify: applies PBC to atoms on a cubic cell"
    print "   -rcut: cutoff distance when building neighbor list, follow by float value"
    print ""

if len(sys.argv) < 2:
    usage()
    exit(0)

#Preprocess Args for flags
rectifyFlag=False
opFlag = False
sliceFlag = False
op = None
rcut = None
lval = 0
for i,v in enumerate(sys.argv):

    if v in ["-rectify","-Rectify"]:
        rectifyFlag = True
        sys.argv.pop(i)

    if v == "-rcut":
        rcut = float(sys.argv[i+1])
        sys.argv.pop(i) #pop the flag
        sys.argv.pop(i) #pop the float

    if v in ["-slice","-Slice"]:
        sliceFlag = True
        sys.argv.pop(i)

#Process args, look for special cases etc.
if sys.argv[1]=="MSD":
    op = "MSD"
    if len(sys.argv)!=5:
        print "Error usage:"
        print "%s MSD POSCAR OUTCAR 25"%sys.argv[0].split("/")[-1]
        exit(0)

    poscar = open(sys.argv[2],"r").readlines()
    outcarFile = sys.argv[3]
    ref = int(sys.argv[4])
    dummy,msd = outcarMeanSquareDisplaceAtom(outcarFile,refStructure=ref)

elif sys.argv[1]=="FILE":
    op = "FILE"
    if len(sys.argv)!=4:
        print "Error usage:"
        print "%s FILE POSCAR OrderParamFile"%sys.argv[0].split("/")[-1]
        exit(0)
    
    poscar = open(sys.argv[2],"r").readlines()
    opfile = sys.argv[3]

    #read in the Order Parameter here, don't parse it later... set flag to false
    opFlag = False
    ops=list()
    for line in open(opfile,"r"):
        if line[0]=="#": continue
        for i in line.split():
            try:
                ops.append(sqrt(float(i)))
            except ValueError:
                pass

elif len(sys.argv)==2:
    poscar=open(sys.argv[1],"r").readlines()
    
elif len(sys.argv)==3:
    op = sys.argv[1]
    if op[:2]=="BO":
        lval=int(op[-1])
        op="BO"
    opFlag = True
    poscar=open(sys.argv[2],"r").readlines()

#Parse POSCAR
[basis,atypes,atoms,head,poscar] = poscarIO.read(poscar)

if rectifyFlag:
    atoms=[[i[0]%basis[0][0],i[1]%basis[1][1],i[2]%basis[2][2]]for i in atoms]

ax,ay,az=zip(*atoms)
if opFlag: ops=[0]*len(ax)
v1,v2,v3=basis
j=0
types=list()
for i in atypes:
    types+=[j+1]*i
    j+=1

fig=mlab.figure()#bgcolor=(1.0,1.0,1.0))

#Get the order parameter and convert to integer format (opsn) for
#coloring of atoms
if opFlag:
    ops,rcut = orderParams[op](np.array(atoms),np.array(basis),l=lval,rcut=rcut)
if op=="MSD":
    ops=np.sqrt(msd.T[-1])

n=None
mnop = min(ops)
mxop = max(ops)
if op=="TET":
    mnop=0.0
    mxop=1.0
    n=11

#Only plot with coloring if there is some variance in the OP
if mxop - mnop > 1E-10:
    #spectral
    mp3d = mlab.points3d(ax,ay,az,ops,colormap='jet',scale_factor=1.9,scale_mode='none',resolution=14)
    if n==None:
        n=min(len(set(ops)),10)
else:
    mp3d = mlab.points3d(ax,ay,az,[mnop]*len(az),colormap='jet',scale_factor=1.9,scale_mode='none',resolution=14)
    n=2

#Color bar formatting
if op != None:
    if op=="MSD":
        op+="^0.5"
    cb = mlab.colorbar(title=op, orientation='vertical', nb_labels=n,nb_colors=n)
    cb.use_default_range = False
    cb.data_range = (mnop,mxop)

#Stupid surrounding box code, sooooo ugly...
z=[0,0,0]
mlab.plot3d([0,v1[0],v1[0]+v2[0],v2[0],0,v3[0]],[0,v1[1],v1[1]+v2[1],v2[1],0,v3[1]],[0,v1[2],v1[2]+v2[2],v2[2],0,v3[2]],color=(0,0,0),line_width=0.5)
mlab.plot3d([v3[0],v3[0]+v1[0],v3[0]+v2[0]+v1[0],v3[0]+v2[0],v3[0]],[v3[1],v3[1]+v1[1],v3[1]+v2[1]+v1[1],v3[1]+v2[1],v3[1]],[v3[2],v3[2]+v1[2],v3[2]+v2[2]+v1[2],v3[2]+v2[2],v3[2]],color=(0,0,0),line_width=0.5)
mlab.plot3d([v1[0],v1[0]+v3[0]],[v1[1],v1[1]+v3[1]],[v1[2],v1[2]+v3[2]],color=(0,0,0),line_width=0.5)
mlab.plot3d([v2[0],v2[0]+v3[0]],[v2[1],v2[1]+v3[1]],[v2[2],v2[2]+v3[2]],color=(0,0,0),line_width=0.5)
mlab.plot3d([v1[0]+v2[0],v1[0]+v2[0]+v3[0]],[v1[1]+v2[1],v1[1]+v2[1]+v3[1]],[v1[2]+v2[2],v1[2]+v2[2]+v3[2]],color=(0,0,0),line_width=0.5)

#average orderParameter based on z-coordinate of atom
if sliceFlag:
    if op==None:
        print "Need an order paramter to plot slice"
        
    #Prepare slices and bins
    zmin = 0
    zmax = basis[2][2]
    nSlices = 21
    delz = (zmax-zmin)/(nSlices-1)
    zBins = [(i)*delz+zmin for i in range(nSlices)]
    zBinLow = [i*delz+zmin for i in range(nSlices-1)]
    zVals = [0.]*nSlices
    zCount = [0]*nSlices

    #Loop over atoms and bin orderParameter based on z-coordinate of atom
    az=[i%zmax for i in az]
    for a,atomZ in enumerate(az):
        i = sum([1 for zLow in zBinLow if atomZ > zLow])
        zCount[i] += 1
        zVals[i] += ops[a]

    #window average and eliminate 0's
    zBins,zVals = map(list,zip(*[[b,z/c] for z,c,b in zip(zVals,zCount,zBins) if c!=0]))
    if zBins[0]!=0:
        
        zBins=[0]+zBins
        zVals=[zVals[-1]]+zVals
    zValsAvg=windowAvg(zVals[-3:]+zVals+zVals[:3],5)[3:-3]

    #plot
    pl.figure()
    pl.plot(zBins,zVals)
    pl.plot(zBins,zValsAvg)
    pl.xlabel("Z")
    pl.ylabel(op)
    pl.xlim(zmin,zmax)
    pl.ylim(0.2,0.8)
    pl.show()

prm.prmshow(fname="POSCAR.png")
