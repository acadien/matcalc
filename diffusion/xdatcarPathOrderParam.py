#!/usr/bin/python

import sys
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy import weave
from scipy.weave import converters
from numpy import *
from mpl_toolkits.mplot3d.art3d import Line3DCollection
#mine
import rootMeanSquareDist
from struct_tools import minImageTranslation

def usage():
    print "%s <xdatfile> <opfile (tetra,cn,rmsd) or \"time\" \"number\"> <comma separated: atom # to render traj>"%(sys.argv[0].split("/")[-1])

if len(sys.argv)!=4:
    usage()
    exit(0)

coarseGrainN=5
xdatcarFileName,opFile,atomsRequested = sys.argv[1],sys.argv[2],map(int,sys.argv[3].split(","))

centralAtomIndex=atomsRequested[0]
atomsRequested=sorted(atomsRequested)
centralAtomIndex=atomsRequested.index(centralAtomIndex)

##Parsing the xdatcar for atomic trajectories
startBlockFlag=False
headlen=0
nAtomsTot = 0
nReq=len(atomsRequested)
atomsTraj = [[] for i in range(nReq)]
reqLookup = dict(zip(atomsRequested,range(nReq)))

for i,line in enumerate(open(xdatcarFileName,"r")):
    if i==0:
        nAtomsTot = int(line.split()[0])
        continue

    if len(line.split())==0:
        headlen=i
        startBlockFlag=True
        continue

    if startBlockFlag and i-headlen-1 in atomsRequested:
        ivd = reqLookup[i-headlen-1]
        atomsTraj[ivd].append(map(float,line.split()))
        if ivd==nReq-1:
            startBlockFlag=False

Ntime = len(atomsTraj[0])
Natom = nReq

##Parsing the orderParam file (can be RMSD, TOP, CN)
jetFlag=False
if opFile in ["time","Time"]:
    jetFlag=True
    atomsOP = [range(Ntime) for i in range(nReq)]
elif opFile in ["number","numb","Number","Numb"]:
    jetFlag=True
    atomsOP = [[float(i)/(nReq-1)]*Ntime for i in range(nReq)]
else:
    atomsOP = [[] for i in range(nReq)]
    for i,line in enumerate(open(opFile,"r")):
        line = line.split()
        if len(line)<20:continue

        for ar,a in reqLookup.iteritems():
            atomsOP[a].append(float(line[ar+1]))

atoms = array(atomsTraj).swapaxes(0,1).ravel()
lengths = array([1.0,1.0,1.0])

weave.inline(rootMeanSquareDist.undoPBCcode,['atoms','Ntime','Natom','lengths'])
atoms.shape = [Ntime,Natom,3]
atomsTraj = atoms.swapaxes(0,1).swapaxes(1,2)

#Apply minimum image distance
centralAtom=atomsTraj[centralAtomIndex,:,0]
for i in range(1,nReq):
    atom=atomsTraj[i,:,0]
    translate = repeat(minImageTranslation(centralAtom,atom,eye(3)).reshape(3,1),atomsTraj.shape[2],axis=1)
    atomsTraj[i,:,:]+=translate

def coarseGrain(x, cg=10):
    return array([sum(x[i*cg:i*cg+cg])/cg for i in range(int(len(x)/cg))]) 

#brg
#winter
#RdBu
#jet
def colorLine(x, y, z, c=None, cmap=plt.get_cmap('cool'),alpha=0.4):
    if c==None:
        c = linspace(0.0, 1.0, len(x))
    c = array(c)
    points = array([x, y, z]).T.reshape(-1, 1, 3)
    segments = concatenate([points[:-1], points[1:]], axis=1)
    lc = Line3DCollection(segments, array=c, cmap=cmap, alpha=alpha)
    return lc
    
fig = plt.figure()
ax = fig.gca(projection='3d')

for atomTraj,atomOP in zip(atomsTraj,atomsOP):
    atomTraj = apply_along_axis(coarseGrain,1,atomTraj,coarseGrainN)
    atomTraj = concatenate((atomTraj[:,0].reshape(3,1),atomTraj,atomTraj[:,-1].reshape(3,1)),axis=1)
    print atomTraj.shape
    if max(atomOP) - min(atomOP) < 1:
        atomOP = concatenate((array([0]),coarseGrain(atomOP,cg=coarseGrainN),array([1])))
    else:
        atomOP = [min(atomOP)] + coarseGrain(atomOP,cg=coarseGrainN) + [max(atomOP)]
    if jetFlag:
        cmap = plt.get_cmap('jet')
    else:
        cmap = plt.get_cmap('cool')
    ax.add_collection(colorLine(atomTraj[0],atomTraj[1],atomTraj[2],c=atomOP,cmap=cmap,alpha=1))
    ax.scatter(atomTraj[0][0],atomTraj[1][0],atomTraj[2][0])

#Set plot bounds so things are scaled poorly
mxx=atomsTraj[:,0,:].max()
mmx=atomsTraj[:,0,:].min()
mxy=atomsTraj[:,1,:].max()
mmy=atomsTraj[:,1,:].min()
mxz=atomsTraj[:,2,:].max()
mmz=atomsTraj[:,2,:].min()
d=max([mxx-mmx,mxy-mmy,mxz-mmz])/2.0
comx=(mxx+mmx)/2
comy=(mxy+mmy)/2
comz=(mxz+mmz)/2

ax.set_xlim(comx-d,comx+d)
ax.set_ylim(comy-d,comy+d)
ax.set_zlim(comz-d,comz+d)

plt.show()
        
