#!/usr/bin/python

import sys
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import meanSquareDist
from scipy import weave
from scipy.weave import converters
from numpy import *

def usage():
    print "%s <xdatfile> <poscarfile> <atom # to render traj>"%(sys.argv[0].split("/")[-1])

if len(sys.argv)!=4:
    usage()
    exit(0)

xdatcarFileName,poscarFile,atomRequested = sys.argv[1],sys.argv[2],int(sys.argv[3])

##Parsing the xdatcar for atomic trajectories
startBlockFlag=False
headlen=0
atomTraj = []
for i,line in enumerate(open(xdatcarFileName,"r")):
    if i==0:
        Nion=int(line.split()[0])
        continue

    if len(line.split())==0:
        headlen=i
        startBlockFlag=True
        continue

    if startBlockFlag and i-headlen-1==atomRequested:
        atomTraj.append(map(float,line.split()))
        startBlockFlag=False

atoms = array([atomTraj]).swapaxes(0,1).ravel()
lengths = array([1.0,1.0,1.0])
Ntime = len(atomTraj)
Natom = 1

weave.inline(meanSquareDist.undoPBCcode,['atoms','Ntime','Natom','lengths'])
atoms.shape = [Ntime,Natom,3]
atomTraj = atoms[:,0,:].swapaxes(0,1)

fig = plt.figure()
ax = fig.gca(projection='3d')

from mpl_toolkits.mplot3d.art3d import Line3DCollection

def coarseGrain(x, y, z,cg=10):
    cgx = [sum(x[i*cg:i*cg+cg])/cg for i in range(int(len(x)/cg))]
    cgy = [sum(y[i*cg:i*cg+cg])/cg for i in range(int(len(x)/cg))]
    cgz = [sum(z[i*cg:i*cg+cg])/cg for i in range(int(len(x)/cg))]
    return cgx,cgy,cgz
atomTraj = coarseGrain(atomTraj[0],atomTraj[1],atomTraj[2],cg=2)

def colorLine(x, y, z, cmap=plt.get_cmap('brg'),alpha=0.9):
    c = linspace(0.0, 1.0, len(x))
    points = array([x, y, z]).T.reshape(-1, 1, 3)
    segments = concatenate([points[:-1], points[1:]], axis=1)
    lc = Line3DCollection(segments, array=c, cmap=cmap, alpha=alpha,lw=3)
    return lc

ax.add_collection(colorLine(atomTraj[0],atomTraj[1],atomTraj[2]))
print min(atomTraj[0])
ax.set_xlim(min(atomTraj[0]),max(atomTraj[0]))
ax.set_ylim(min(atomTraj[1]),max(atomTraj[1]))
ax.set_zlim(min(atomTraj[2]),max(atomTraj[2]))
#mn,mx = min(map(min,[atomTraj[0],atomTraj[1],atomTraj[2]])), max(map(max,[atomTraj[0],atomTraj[1],atomTraj[2]]))
#ax.set_xlim(mn,mx)
#ax.set_ylim(mn,mx)
#ax.set_zlim(mn,mx)
xt,yt,zt = ax.get_xticks(),ax.get_yticks(),ax.get_zticks()
ax.set_xticklabels(["" for i in xt])
ax.set_yticklabels(["" for i in yt])
ax.set_zticklabels(["" for i in zt])
plt.show()
        
