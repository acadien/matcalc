#!/usr/bin/python

import sys
from numpy import *
#mine
from datatools import windowAvg

def usage():
    print "%s <hopped.dat> <outcar.rmsd> <opt:opfile (tetra,cn,rmsd)>"%(sys.argv[0].split("/")[-1])
    
if len(sys.argv) not in [3,4]:
    usage()
    exit(0)

parseHop= lambda x: map(float,x.split()) 
hoppedfile = sys.argv[1]
try:
    atom,jump50,jump100,jump200,jump500,jump1000 = zip(*map(parseHop,open(hoppedfile,"r").readlines()))
except ValueError:
    atom,jump50,jump100,jump200,jump500,jump1000 = zip(*map(parseHop,open(hoppedfile,"r").readlines()[1:]))
atomsJumped200 = array([i for i,j in enumerate(jump200) if j>0])

rmsdFile = sys.argv[2]
if rmsdFile[-4:]!="rmsd": exit(0)

#use jump50, 50ts for 2.5A seems to be a good param for capturing hops/flow.
rmsdAvg=list()
rmsdPerAtom=list()
for i,line in enumerate(open(rmsdFile,"r").readlines()):
    if i==0:
        continue
    line = map(lambda x: float(x),line.split())
    rmsdAvg.append(line[0])
    rmsdPerAtom.append(array(line[1:])[atomsJumped200])
rmsdPerAtom=array(rmsdPerAtom)

##Parsing the orderParam file (can be RMSD, TOP, CN)
atomOP = []
if len(sys.argv)==4:
    opFile=sys.argv[3]
    for i,line in enumerate(open(opFile,"r")):
        line = line.split()
        if len(line)<20:
            continue
        atomOP.append(array(map(float,line[1:]))[atomsJumped200])
atomOP=array(atomOP)

import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import Grid
from matplotlib import collections as mc
fig=plt.figure()

mxx = rmsdPerAtom.max()*1.1

#figure out the geometry of the plotting window
natom = rmsdPerAtom.shape[1]
nr,nc = 25,25
print natom
if natom<49:
    nr,nc = 6,8
elif natom<97:
    nr,nc = 12,8
elif natom<105:
    nr,nc = 13,8
elif natom<151:
    nr,nc = 15,10
elif natom<193:
    nr,nc = 16,12
elif natom<253:
    nr,nc = 18,14
elif natom<289:
    nr,nc = 18,16
grid = Grid(fig, rect=111, nrows_ncols=(nr,nc),
            axes_pad=0.05, label_mode='L')

def coarseGrain(x, cg=10):
    return [sum(x[i*cg:i*cg+cg])/cg for i in range(int(len(x)/cg))]

def colorLine(x, y, c=None, cmap=plt.get_cmap('cool'),alpha=0.4):
    if c==None:
        c = linspace(0.0, 1.0, len(x))
    c = array(c)
    points = array([x, y]).T.reshape(-1, 1, 2)
    segments = concatenate([points[:-1], points[1:]], axis=1)
    lc = mc.LineCollection(segments, array=c, cmap=cmap, alpha=alpha)
    return lc

atomOP=swapaxes(atomOP,0,1)
for i,rmsd in enumerate(rmsdPerAtom.swapaxes(0,1)):
    if len(atomOP)==0:            
        grid[i].plot(coarseGrain(rmsd,100),label = str(atomsJumped200[i]))
    else:
        cgy,cgv = coarseGrain(rmsd,100),coarseGrain(atomOP[i],100)
        cgx=linspace(0.0,len(rmsd)-1,len(cgy))
        grid[i].add_collection(colorLine(cgx,cgy,c=cgv))
    grid[i].set_xlim([0,cgx[-1]])
    grid[i].set_ylim([0,mxx])
    grid[i].text(0.5,0.8,str(atomsJumped200[i]),horizontalalignment='center',verticalalignment='center',transform=grid[i].transAxes)
    
plt.tight_layout()
plt.show()
