#!/usr/bin/python

import plotRemoteMaya as prm #mine

import sys
import pylab as pl
from numpy import array
from mayavi import mlab

def usage():
    print "%s <dump file with atomic coordinates> <start config>"%sys.argv[0].split("/")[-1]

if len(sys.argv)<2:
    usage()
    exit(0)

#Prepare data
fname=sys.argv[1]
fraw = open(fname,"r").readlines()
configCount=0
for i,line in enumerate(fraw):
    if "ITEM: TIMESTEP" in line:
        configCount +=1

startConfig=0
if len(sys.argv)>=3:
    startConfig=int(sys.argv[2])

#Prepare plotting
mxc=20
iline=0
natoms=0
plotForces=False
c=0
while True:
    cfgFound=False
    info="Configuration Number %d\n"%c
    for i,line in enumerate(fraw[iline:]):
        if "ITEM: TIMESTEP" in line:
            info+="Timestep = %d\n"%int(fraw[iline+i+1])
        if "ITEM: NUMBER OF ATOMS" in line:
            natoms=int(fraw[iline+i+1])
            info+="N-Atoms = %d\n"%natoms 
        if "ITEM: ATOMS" in line:
            if "fx" in line:
                plotForces=True
            cfgFound=True
            iline+=i+1
            break
    
    #no more configurations, stop plotting yo
    if not cfgFound: break
    c+=1
    if c < startConfig:
        continue


    print c,"/",configCount

    if not plotForces:
        #Parsing function for parsing atomic coordinates and types
        parse=lambda x:[int(x[0]),float(x[1]),float(x[2]),float(x[3]),float(x[4])]
        #Apply the parse on the portion of the datafile containing coordinates
        number,colors,xs,ys,zs = \
            zip(*[parse(line.split()) for line in fraw[iline:iline+natoms]])
    else:
        #Parsing function for parsing atomic coordinates and types
        parse=lambda x:[int(x[0]),float(x[1]),float(x[2]),float(x[3]),float(x[4]),float(x[5]),float(x[6]),float(x[7])]
        #Apply the parse on the portion of the datafile containing coordinates
        number,colors,xs,ys,zs,fx,fy,fz = \
            zip(*[parse(line.split()) for line in fraw[iline:iline+natoms]])

    iline+=natoms
    
    cmn=min(colors)
    cmx=max(colors)

    xmn,xmx=min(xs),max(xs)
    ymn,ymx=min(ys),max(ys)
    zmn,zmx=min(zs),max(zs)

    xs = array(xs)
    ys = array(ys)
    zs = array(zs)

    if (xs.sum()+ys.sum()+zs.sum())/3/len(xs) < 1.0:
        xs*=(xmx-xmn)
        ys*=(ymx-ymn)
        zs*=(zmx-zmn)

    fig=mlab.figure(bgcolor=(0.0,0.0,0.0))

    if plotForces:
        mlab.quiver3d(xs,ys,zs,fx,fy,fz,scale_factor=1.0,scale_mode='none')
        mlab.colorbar(title="intensity")
    else:
        mlab.points3d(xs,ys,zs,[1.0]*len(zs),scale_factor=0.1,scale_mode='none')
        
    mlab.show()
#prm.prmshow(fname="blasni.png")
