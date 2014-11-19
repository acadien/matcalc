#!/usr/bin/python

import sys
#mine
import plotRemoteMaya as prm
import plotRemote as pr

from mayavi import mlab
import numpy as np
import pylab as pl
from math import *
#mine
from orderParam import coordinationNumber,bondOrientation,tetrahedral
from outcarPlotMSDAtom import outcarMeanSquareDisplaceAtom
import datatools
from neighbors import neighbors
import poscarIO,lammpsIO,outcarIO

orderParams={"CN":coordinationNumber, \
             "BO":bondOrientation, \
             "TET":tetrahedral, \
             "TET2":tetrahedral}

#RMSD is not an order parameter but useful for labeling atoms at high temperatures.
#FILE is not an order parameter but applies colors based data from a file

def usage():
    print "----------------------------------------------------------------------------"
    print "%s <order parameter> <POSCAR/dump file> <OPT: OUTCAR, REF>"%sys.argv[0]
    print "----------------------------------------------------------------------------"
    print "Order Parameter can be one of:"
    print "     CN : Coordination Number"
    print "    BO# : Bond Orientation (Q) with l=#"
    print "    TET : Tetrahedral order parameter (Sg)"
    print "   TET2 : 2nd shell averaged TET (Sg)"
    print "   RMSD : Mean Square Displacement, requires OUTCAR (annealed MD simulation)"
    print "   FILE : Parses a file (column, row, or CSV data) and applies coloring from that data"
    print "----------------------------------------------------------------------------"
    print "Flags can be used anywhere in args:"
    print "-rectify  : applies PBC to atoms on a cubic cell"
    print "   -rcut #: cutoff distance when building neighbor list, follow by float value"
    print "   -hist  : generates a histogram of the order parameter"
    print "   -N #   : selects a configuration to use"
    print "-bounds #,# : min,max bounds on order parameter"
    print ""

if len(sys.argv) < 2:
    usage()
    exit(0)

#Preprocess Args for flags
rectifyFlag=False
opFlag = False
sliceFlag = False
histFlag = False
op = None
rcut = None
Nconfig = 0
lval = 0
toPop=list()
minv,maxv = None,None
for i,v in enumerate(sys.argv):

    if v in ["-rectify","-Rectify"]:
        rectifyFlag = True
        toPop.append(i)

    if v == "-rcut":
        rcut = float(sys.argv[i+1])
        toPop.append(i)
        toPop.append(i+1)

    if v in ["-slice","-Slice"]:
        sliceFlag = True
        toPop.append(i)

    if v in ["-hist","-Hist","-HIST"]:
        histFlag = True
        toPop.append(i)

    if v in ["-N","-n"]:
        Nconfig = int(sys.argv[i+1])
        toPop.append(i)
        toPop.append(i+1)

    if v in ["-h"]:
        usage()
        exit(0)

    if v in ["-bound","-bounds","-bnds","-bnd"]:
        minv,maxv = map(float,sys.argv[i+1].split(","))
        toPop.append(i)
        toPop.append(i+1)

sys.argv = [sys.argv[i] for i in range(len(sys.argv)) if i not in toPop]

#Process args, look for special cases etc.
if sys.argv[1]=="RMSD":
    op = "RMSD"
    if len(sys.argv)!=5:
        print "Error usage:"
        print "%s RMSD POSCAR OUTCAR 25"%sys.argv[0].split("/")[-1]
        exit(0)

    configFile = sys.argv[2]
    poscar = open(configFile,"r").readlines()
    outcarFile = sys.argv[3]
    ref = int(sys.argv[4])
    dummy,msd = outcarMeanSquareDisplaceAtom(outcarFile,refStructure=ref)

elif sys.argv[1]=="FILE":
    op = "FILE"
    if len(sys.argv)!=4:
        print "Error usage:"
        print "%s FILE POSCAR OrderParamFile"%sys.argv[0].split("/")[-1]
        exit(0)
    
    configFile = sys.argv[2]
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

elif len(sys.argv) in [2,4]:
    configFile=sys.argv[1]
    
elif len(sys.argv)==3:
    op = sys.argv[1]
    if op[:2]=="BO":
        lval=int(op[-1])
        op="BO"
    opFlag = True
    configFile=sys.argv[2]

else:
    usage()
    exit(0)

#Parse POSCAR
j=0
if "POSCAR" in configFile or "CONTCAR" in configFile:
    poscar=open(configFile,"r").readlines()
    [basis,atypes,atoms,head,poscar] = poscarIO.read(poscar)
    types=list()
    for i in atypes:
        types+=[j+1]*i
        j+=1

#Parse OUTCAR
elif "OUTCAR" in configFile:
    TE,stress,basis,atoms,forces,types = outcarIO.outcarReadConfig(configFile,Nconfig)

#Parse LAMMPS
else: 
    basis,types,atoms = lammpsIO.readConfig(configFile,Nconfig)

ax,ay,az=map(np.array,zip(*atoms))
nAtoms = len(ax)
v1,v2,v3=basis


if rectifyFlag:
    ax,ay,az = ax%v1[0], ay%v2[1], az%v3[2]
    atoms = np.vstack((ax,ay,az))
    atoms.shape = nAtoms,3

#Set default rcut value for tetrahedral ordering
if op in ["TET","TET2"] and rcut==None:
    rcut=3.1

#Get the order parameter and convert to integer format (opsn) for coloring of atoms
if opFlag:
    ops,rcut = orderParams[op](np.array(atoms),np.array(basis),l=lval,rcut=rcut)
else:
    ops = types


if op=="RMSD":
    ops = np.sqrt(msd.T[-1])

if op=="TET2": 
    bounds = [[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
    neighbs = neighbors(atoms,bounds,rcut)
    secondShell = list()

    for a,firstNeighbs in enumerate(neighbs):
        secondNeighbs = [neighbs[i] for i in firstNeighbs]
        totalNeighb = set(datatools.flatten([firstNeighbs]+secondNeighbs))
        N = len(totalNeighb)

        if N==0:
            secondShell.append(ops[a])
        else:
            secondShell.append(sum([ops[i] for i in totalNeighb])/N)

    ops = secondShell

n=None
mnop = min(ops)
mxop = max(ops)
if op in ["TET","TET2"]:
    mnop, mxop, n = 0.0, 1.0, 11

#Auto-set the resolution so you don't burn your computer down trying to render this
nAtom = len(ax)
res=16.0
while nAtom > 1100 and res>3.0:
    nAtom/=4.0
    res/=2.0
res=max(res,3.0)
nAtom = len(ax)

if histFlag:
    if op==None:
        print "Need an order parameter in order to generate histogram"
        exit(0)

    pl.figure()
    pl.hist(ops)
    pl.xlabel(op)
    pl.ylabel("Count")
    pl.draw() #needed when using interactive plotting
    pr.prshow()
    exit(0)

fig=mlab.figure(bgcolor=(0.8,0.8,0.8))

#Only plot with coloring if there is some variance in the OP
if mxop - mnop > 1E-10:
    #spectral
    mp3d = mlab.points3d(ax,ay,az,ops,colormap='jet',scale_factor=1.9,scale_mode='none',resolution=res)
    if n==None:
        n=min(len(set(ops)),10)
else:
    mp3d = mlab.points3d(ax,ay,az,[mnop]*len(az),colormap='jet',scale_factor=1.9,scale_mode='none',resolution=res)
    n=2

#Color bar formatting
if op != None:
    if op=="RMSD":
        op+="^0.5"
    cb = mlab.colorbar(title=op, orientation='vertical', nb_labels=n,nb_colors=n)
    cb.use_default_range = False

    if minv == None and maxv==None:
        cb.data_range = (mnop,mxop)
    else:
        cb.data_range = (minv,maxv)


#Stupid surrounding box code, sooooo ugly...
z=[0,0,0]
mlab.plot3d([0,v1[0],v1[0]+v2[0],v2[0],0,v3[0]],[0,v1[1],v1[1]+v2[1],v2[1],0,v3[1]],[0,v1[2],v1[2]+v2[2],v2[2],0,v3[2]],color=(0,0,0),line_width=0.5)
mlab.plot3d([v3[0],v3[0]+v1[0],v3[0]+v2[0]+v1[0],v3[0]+v2[0],v3[0]],[v3[1],v3[1]+v1[1],v3[1]+v2[1]+v1[1],v3[1]+v2[1],v3[1]],[v3[2],v3[2]+v1[2],v3[2]+v2[2]+v1[2],v3[2]+v2[2],v3[2]],color=(0,0,0),line_width=0.5)
mlab.plot3d([v1[0],v1[0]+v3[0]],[v1[1],v1[1]+v3[1]],[v1[2],v1[2]+v3[2]],color=(0,0,0),line_width=0.5)
mlab.plot3d([v2[0],v2[0]+v3[0]],[v2[1],v2[1]+v3[1]],[v2[2],v2[2]+v3[2]],color=(0,0,0),line_width=0.5)
mlab.plot3d([v1[0]+v2[0],v1[0]+v2[0]+v3[0]],[v1[1]+v2[1],v1[1]+v2[1]+v3[1]],[v1[2]+v2[2],v1[2]+v2[2]+v3[2]],color=(0,0,0),line_width=0.5)

pl.ion() #turn on interactive plotting

#average orderParameter based on z-coordinate of atom
if sliceFlag:
    if op==None:
        print "Need an order parameter to plot slice"
        exit(0)

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
    zValsAvg=datatools.windowAvg(zVals[-3:]+zVals+zVals[:3],5)[3:-3]

    #plot
    pl.figure()
    pl.plot(zBins,zVals)
    pl.plot(zBins,zValsAvg)
    pl.xlabel("Z")
    pl.ylabel(op)
    pl.xlim(zmin,zmax)
    pl.draw()
    pl.show()

prm.prmshow(fname="plot3.png")
