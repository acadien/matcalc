#!/usr/bin/python

import sys
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import meanSquareDist
from scipy import weave
from scipy.weave import converters
from numpy import *
from mpl_toolkits.mplot3d import proj3d
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

fig = plt.figure(figsize=(3,3))
ax = fig.gca(projection='3d')
ax.view_init(elev=-45., azim=-45)
from mpl_toolkits.mplot3d.art3d import Line3DCollection

def coarseGrain(x, y, z,cg=4):
    cgx = [sum(x[i*cg:i*cg+cg])/cg for i in range(int(len(x)/cg))]
    cgy = [sum(y[i*cg:i*cg+cg])/cg for i in range(int(len(x)/cg))]
    cgz = [sum(z[i*cg:i*cg+cg])/cg for i in range(int(len(x)/cg))]
    return cgx,cgy,cgz
cut=-1
atomTraj = coarseGrain(atomTraj[0][:cut],atomTraj[1][:cut],atomTraj[2][:cut],cg=1)

nx,xx = min(atomTraj[0]),max(atomTraj[0])
ny,xy = min(atomTraj[1]),max(atomTraj[1])
nz,xz = min(atomTraj[2]),max(atomTraj[2])

print xz-nz

dx=max(((xx-nx)/2.0,(xy-ny)/2.0,(xz-nz)/2.0))

cx=(xx+nx)/2.0
cy=(xy+ny)/2.0
cz=(xz+nz)/2.0

nx,xx=cx-dx,cx+dx
ny,xy=cy-dx,cy+dx
nz,xz=cz-dx,cz+dx

#nx,xx = min(atomTraj[0]),max(atomTraj[0])
#ny,xy = min(atomTraj[1]),max(atomTraj[1])
#nz,xz = min(atomTraj[2]),max(atomTraj[2])

ax.set_xlim(nx,xx)
ax.set_ylim(ny,xy)
ax.set_zlim(nz,xz)
#ax.set_zlim(nz,xz)


ax.plot([nx,xx],[ny,ny],zs=[nz,nz],color="white",lw=1,zorder=9)
ax.plot([nx,xx],[xy,xy],zs=[nz,nz],color="white",lw=1)
ax.plot([nx,xx],[ny,ny],zs=[xz,xz],color="white",lw=1)
ax.plot([nx,xx],[xy,xy],zs=[xz,xz],color="white",lw=1)

ax.plot([nx,nx],[ny,xy],zs=[nz,nz],color="white",lw=1)
ax.plot([xx,xx],[ny,xy],zs=[nz,nz],color="white",lw=1,zorder=0)
ax.plot([nx,nx],[ny,xy],zs=[xz,xz],color="white",lw=1)
ax.plot([xx,xx],[ny,xy],zs=[xz,xz],color="white",lw=1)

ax.plot([nx,nx],[ny,ny],zs=[nz,xz],color="white",lw=1)
ax.plot([xx,xx],[ny,ny],zs=[nz,xz],color="white",lw=1)
ax.plot([nx,nx],[xy,xy],zs=[nz,xz],color="white",lw=1)
ax.plot([xx,xx],[xy,xy],zs=[nz,xz],color="white",lw=1)


def colorLine(x, y, z, cmap=plt.get_cmap('gnuplot'),alpha=1.0):
    c = linspace(0.0, 1.0, len(x))
    points = array([x, y, z]).T.reshape(-1, 1, 3)
    segments = concatenate([points[:-1], points[1:]], axis=1)
    lc = Line3DCollection(segments, array=c, cmap=cmap, alpha=alpha,lw=2)
    return lc

ax.add_collection(colorLine(atomTraj[0],atomTraj[1],atomTraj[2]))
#ax.w_xaxis.set_pane_color((0.0,0.0,0.0,1.0))
#ax.w_yaxis.set_pane_color((0.0,0.0,0.0,1.0))
#ax.w_zaxis.set_pane_color((0.0,0.0,0.0,1.0))
ax.w_xaxis.set_pane_color((0.1,0.1,0.1,1.0))
ax.w_yaxis.set_pane_color((0.1,0.1,0.1,1.0))
ax.w_zaxis.set_pane_color((0.1,0.1,0.1,1.0))

#n,mx = min(map(min,[atomTraj[0],atomTraj[1],atomTraj[2]])), max(map(max,[atomTraj[0],atomTraj[1],atomTraj[2]]))
#ax.set_xlim(n,x)
#ax.set_ylim(n,x)
#ax.set_zlim(n,x)
xt,yt,zt = ax.get_xticks(),ax.get_yticks(),ax.get_zticks()
ax.set_xticks([])#labels(["" for i in xt])
ax.set_yticks([])#labels(["" for i in yt])
ax.set_zticks([])#labels(["" for i in zt])


def orthogonal_proj(zfront, zback):
    a = (zfront+zback)/(zfront-zback)
    b = -2*(zfront*zback)/(zfront-zback)
    # -0.0001 added for numerical stability as suggested in:
    # http://stackoverflow.com/questions/23840756
    return array([[1,0,0,0],
                        [0,1,0,0],
                        [0,0,a,b],
                        [0,0,-0.0001,zback]])
 
# Later in your plotting code ...
proj3d.persp_transformation = orthogonal_proj

#plt.savefig("/home/acadien/Desktop/liq.svg") #63 XDATCAR CONTCAR
#plt.savefig("/home/acadien/Desktop/hop.svg") #65 XDATCAR0 CONTCAR0
plt.show()
        
