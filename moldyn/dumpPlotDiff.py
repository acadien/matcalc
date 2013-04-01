#!/usr/bin/python

import plotRemoteMaya as prm #mine

import argparse
import pylab as pl
from numpy import array
from mayavi import mlab

parser = argparse.ArgumentParser(description="Plots atoms distance between them at differet iterations from a single dump file, configurations are 1 indexed")
parser.add_argument('filename',type=str,nargs=1)
parser.add_argument('configA',type=int,nargs=1)
parser.add_argument('configB',type=int,nargs=1)

args = parser.parse_args()
fname = args.filename[0]
configA = args.configA[0]
configB = args.configB[0]

fraw = open(fname,"r").readlines()
configCount=0
for i,line in enumerate(fraw):
    if "ITEM: TIMESTEP" in line:
        configCount +=1

if configA == 0: configA = 1
if configB == -1: configB = configCount

if configA > configCount or configB > configCount:
    print "Invalid configuration selection, max config = %d."%configCount
    exit(0)

#Prepare plotting
mxc=20
iline=0
natoms=0
c=0
atomsA=list()
atomsB=list()
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
            cfgFound=True
            iline+=i+1
            break
    
    #no more configurations, stop plotting yo
    if not cfgFound: break
    
    c+=1
    if c not in [configA,configB]:
        continue

    if c == configA:
        parse=lambda x:[int(x[0]),float(x[1]),float(x[2]),float(x[3]),float(x[4])]
        #Apply the parse on the portion of the datafile containing coordinates
        number,colors,xs,ys,zs = \
            zip(*[parse(line.split()) for line in fraw[iline:iline+natoms]])

        iline+=natoms
    
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

        atomsA = [xs,ys,zs]

    if c == configB:
        parse=lambda x:[int(x[0]),float(x[1]),float(x[2]),float(x[3]),float(x[4])]
        #Apply the parse on the portion of the datafile containing coordinates
        number,colors,xs,ys,zs = \
            zip(*[parse(line.split()) for line in fraw[iline:iline+natoms]])

        iline+=natoms
    
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

        atomsB = [xs,ys,zs]

fig=mlab.figure(bgcolor=(1.0,1.0,1.0))

[xa,ya,za]=atomsA
[xb,yb,zb]=atomsB

d=((xb-xa)**2+(yb-ya)**2+(zb-za)**2)**0.5/3

mlab.quiver3d(xa,ya,za,xb-xa,yb-ya,zb-za,scalars=d,scale_factor=1.0,scale_mode='none',line_width=4.0)
cb = mlab.colorbar(title="Distance in Angstroms",)
cb._title_text_property.color = (0.0,0.0,0.0)
cb._label_text_property.color = (0.0,0.0,0.0)
mlab.show()

