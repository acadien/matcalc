#!/usr/bin/python

#mine
import plotRemote
#theirs
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
import utils

utils.usage(["<xdatfile>","opfile (tetra,cn,rmsd or \"time\" or \"number\""],2,2)

coarseGrainN=1
xdatcarFileName,opFile = sys.argv[1],sys.argv[2]

#centralAtomIndex=atomsRequested[0]
#centralAtomIndex=atomsRequested.index(centralAtomIndex)

##Parsing the xdatcar for atomic trajectories
startBlockFlag=False
headlen=0
nAtom = 0
atomsTraj = []
count = 0
for i,line in enumerate(open(xdatcarFileName,"r")):
    if i==0:
        nAtom = int(line.split()[0])
        atomsTraj = [[] for i in range(nAtom)]
        continue

    if len(line.split())==0:
        headlen=i
        startBlockFlag=True
        continue

    if startBlockFlag:
        ivd = i-headlen-1
        atomsTraj[ivd].append(map(lambda x:(float(x)+0.5)%1.0,line.split()))
        if ivd==nAtom-1:
            count+=1
            startBlockFlag=False
    if count==1000:
        break

nTime = len(atomsTraj[0])
print nTime
nAtom = nAtom

##Parsing the orderParam file (can be RMSD, TOP, CN)
jetFlag=False
if opFile in ["time","Time"]:
    jetFlag=True
    atomsOP = [range(nTime) for i in range(nAtom)]
elif opFile in ["number","numb","Number","Numb"]:
    jetFlag=True
    atomsOP = [[float(i)/(nAtom-1)]*nTime for i in range(nAtom)]
else:
    atomsOP = []#[[] for i in range(nAtom)]
    count = 0
    for i,line in enumerate(open(opFile,"r")):
        line = line.split()
        if len(line)<20:continue

        count+=1
        atomsOP.append(map(float,line[1:]))
        if count==1000:
            break
    atomsOP = zip(*atomsOP)
atoms = array(atomsTraj).swapaxes(0,1).ravel()
lengths = array([1.0,1.0,1.0])
b = eye(3).ravel()*19.058821406

weave.inline(rootMeanSquareDist.undoPBCcode,['atoms','nTime','nAtom','b'])
atoms.shape = [nTime,nAtom,3]
atomsTraj = atoms.swapaxes(0,1).swapaxes(1,2)

#Apply minimum image distance
#centralAtom=atomsTraj[centralAtomIndex,:,0]
#for i in range(1,nAtom):
#    atom=atomsTraj[i,:,0]
#    translate = repeat(minImageTranslation(centralAtom,atom,eye(3)).reshape(3,1),atomsTraj.shape[2],axis=1)
#    atomsTraj[i,:,:]+=translate

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
    print x[:10]
    print y[:10]
    print z[:10]
    points = array([x, y, z]).T.reshape(-1, 1, 3)
    segments = concatenate([points[:-1], points[1:]], axis=1)
    alphas = linspace(0.0, 1.0, len(x))**2
    rgba = cmap(c)
    rgba[:,3] = alphas
    lc = Line3DCollection(segments, edgecolor = rgba)
    return lc
    
fig = plt.figure(figsize=(18,18))
ax = fig.gca(projection='3d')

for atomTraj,atomOP in zip(atomsTraj,atomsOP):
    #atomTraj = apply_along_axis(coarseGrain,1,atomTraj,coarseGrainN)
    #atomTraj = concatenate((atomTraj[:,0].reshape(3,1),atomTraj,atomTraj[:,-1].reshape(3,1)),axis=1)
    #if max(atomOP) - min(atomOP) < 1:
    #    atomOP = concatenate((array([0]),coarseGrain(atomOP,cg=coarseGrainN),array([1])))
    #else:
    #    atomOP = [min(atomOP)] + coarseGrain(atomOP,cg=coarseGrainN) + [max(atomOP)]
    if jetFlag:
        cmap = plt.get_cmap('jet')
    else:
        cmap = plt.get_cmap('cool')
    if min(atomOP)<0.7:
        ax.add_collection(colorLine(atomTraj[0],atomTraj[1],atomTraj[2],c=atomOP,cmap=cmap,alpha=1))
        break
    #ax.scatter(atomTraj[0][0],atomTraj[1][0],atomTraj[2][0])

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
#ax.savefig("temp.png")
plotRemote.prshow()
        
