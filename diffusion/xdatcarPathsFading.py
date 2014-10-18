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

nTime = len(atomsTraj[0])
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

    atomsOP = zip(*atomsOP)

atoms = array(atomsTraj).swapaxes(0,1)
basis = eye(3)

atoms = rootMeanSquareDist.unwrap(atoms,basis)
atomsTraj = atoms.swapaxes(0,1).swapaxes(1,2)

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
    alphas = linspace(0.0, 1.0, len(x))**2
    cmap = plt.get_cmap('Blues')
    rgba = cmap(c)
    rgba[:,3] = alphas
    lw = linspace(16,3,len(x))
    lc = Line3DCollection(segments, edgecolor = rgba, lw=lw)
    return lc
    
fig = plt.figure(figsize=(18,18))

mxx=atomsTraj[:,0,0].max()
mmx=atomsTraj[:,0,0].min()
mxy=atomsTraj[:,1,0].max()
mmy=atomsTraj[:,1,0].min()
mxz=atomsTraj[:,2,0].max()
mmz=atomsTraj[:,2,0].min()
d=max([mxx-mmx,mxy-mmy,mxz-mmz])/2.0
comx=(mxx+mmx)/2
comy=(mxy+mmy)/2
comz=(mxz+mmz)/2

ax = fig.gca(projection='3d',frame_on=False)
ax.set_frame_on(False)
#ax.set_axes([comx-d,comy-d,comz-d,comx+d,comy+d,comy+d])
plt.tight_layout(pad=-1.0,h_pad=-1.0,w_pad=-1.0)

l=100
m=0
for a in range(0,nTime-l,5):
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
        if max(atomOP)>1.0:
            ax.add_collection(colorLine(atomTraj[0][a:a+l],atomTraj[1][a:a+l],atomTraj[2][a:a+l],c=atomOP[a:a+l],cmap=cmap,alpha=1))
        break
    ax.set_frame_on(False)
    ax.set_xticklabels([""]*len(ax.get_xticklabels()),visible=False)
    ax.set_yticklabels([""]*len(ax.get_yticklabels()),visible=False)
    ax.set_zticklabels([""]*len(ax.get_zticklabels()),visible=False)
    ax.set_xlim(atomTraj[0].min(),atomTraj[0].max())
    ax.set_ylim(atomTraj[1].min(),atomTraj[1].max())
    ax.set_zlim(atomTraj[2].min(),atomTraj[2].max())
#        ax.set_xlim(comx-d,comx+d)
#        ax.set_ylim(comy-d,comy+d)
#        ax.set_zlim(comz-d,comz+d)

    ax.set_frame_on(False)
    fig.savefig("%4.4d.png"%m)#,bbox_inches='tight',pad_inches=0.0,edgecolor="red",facecolor="blue")
    m+=1
    plt.cla()
    print a

    #plotRemote.prshow()
        
