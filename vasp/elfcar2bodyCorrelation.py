#!/usr/bin/python

import chgcarIO
import sys
import pylab as pl
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
import itertools

def usage():
    print "%s <charge density> <rval>"%sys.argv[0].split("/")[-1]
    print "%s <charge density> <rstart> <rstop> <rsteps>"%sys.argv[0].split("/")[-1]
    print "calculates 2 body spatial correlation of a charge density file"

#Returns an array of grid points that lay within a sphere of radius rCut, grid is on rDel
def sphereCoords(delR,rCut):
    rStep=int(np.ceil(rCut/delR))
    rCut2 =rCut*rCut

    #inefficient loop who cares rStep is small, this func is used only once
    coords=list()
    for i in range(-1*rStep,rStep+1):
        for j in range(-1*rStep,rStep+1):
            for k in range(-1*rStep,rStep+1):
                if i*i+j*j+k*k < rCut2:
                    coords.append(np.array([i,j,k]))
    return np.array(coords)

#given a central point, shift coordinates and boundary, get the PBC coordinates
def gen2BCoords(center,coordShifts,bounds):
    return np.mod( center + coordShifts, bounds)

if len(sys.argv)<2:
    usage()
    exit(0)

rVals=None
if len(sys.argv)==3:
    chgcarfile = sys.argv[1]
    rVals = [float(sys.argv[2])]
if len(sys.argv)==5:
    chgcarfile = sys.argv[1]
    rStart = float(sys.argv[2])
    rStop = float(sys.argv[3])
    rSteps = int(sys.argv[4])
    rVals = [(rStop-rStart)*i/rSteps+rStart for i in range(rSteps)]

chgcar=open(chgcarfile,"r").readlines()
pcar,chg = chgcarIO.read(chgcar)
(basis,atypes,atoms,head)=pcar

#R space distance between grid points
delR = sum([np.linalg.norm(basis[i])/chg.shape[i] for i in range(3)])/3.0

bounds=chg.shape
xGrids,yGrids,zGrids=[range(i) for i in bounds]
gridCoords=[i for i in itertools.product(xGrids,yGrids,zGrids)]
for rCut in rVals:
    chgs=list()
    coordShifts=sphereCoords(delR,rCut)
    N=coordShifts.shape[0]
    a=0
    for gc in gridCoords:
        x,y,z=np.split(gen2BCoords(gc,coordShifts,bounds),3,axis=1)
        chgs.append(chg[x,y,z].sum()/N)
        a+=1
        if a%10000==0:
            print a
    vals,bins=np.histogram(chgs,100)
    datatowrite=["chg_density count\n"]+map(lambda x:" ".join(map(str,x))+"\n",zip(bins,vals))
    open(chgcarfile+"2body.RCut%s"%str(round(rCut,2)),"w").writelines(datatowrite)
