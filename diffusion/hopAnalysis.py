#!/usr/bin/python

import sys
from numpy import *

parseHop = lambda x: map(float,x.split())
atom,jump50,jump100,jump200 = zip(*map(parseHop,open("hoppedAtomsOCFULL.dat","r").readlines()))
atomsJumped200 = array([i for i,j in enumerate(jump200) if j>0])

rmsdFile = sys.argv[1]
if rmsdFile[-4:]!="rmsd": exit(0)

#use jump50, 50ts for 2.5A seems to be a good param for capturing hops/flow.
rmsdAvg=list()
rmsdPerAtom=list()
for i,line in enumerate(open(rmsdFile,"r").readlines()):
    if i==0:
        continue
    line = map(lambda x: sqrt(float(x)),line.split())
    rmsdAvg.append(line[0])
    rmsdPerAtom.append(array(line[1:])[atomsJumped200])
rmsdPerAtom=array(rmsdPerAtom)

import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import Grid
fig=plt.figure()
grid = Grid(fig, rect=111, nrows_ncols=(11,6),
            axes_pad=0.05, label_mode='L')
#fig, axesGrid = plt.subplots(nrows=8,ncols=11)
#axesList = [axes for subaxes in axesGrid for axes in subaxes]
for i,v in enumerate(rmsdPerAtom.swapaxes(0,1)):
    grid[i].plot(v,label = str(atomsJumped200[i]))
    grid[i].set_xticklabels([])
    grid[i].set_ylim([0,3.5])
    grid[i].text(0.9,0.9,str(atomsJumped200[i]),horizontalalignment='center',verticalalignment='center',transform=grid[i].transAxes)
plt.show()
