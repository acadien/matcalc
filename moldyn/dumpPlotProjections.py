#!/usr/bin/python

import plotRemote as pr #mine

import sys
import pylab as pl

def usage():
    print "%s <dump file with atomic coordinates>"%sys.argv[0].split("/")[-1]

if len(sys.argv)!=2:
    usage()
    exit(0)

#Prepare data
fname=sys.argv[1]
fraw = open(fname,"r").readlines()

#Prepare plotting

mxc=20
iline=0
natoms=0
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
            cfgFound=True
            iline+=i+1
            break
    
    #no more configurations, stop plotting yo
    if not cfgFound: break
    c+=1

    #Parsing function for parsing atomic coordinates and types
    parse=lambda x:[int(x[0]),float(x[1])/mxc,float(x[2]),float(x[3]),float(x[4])]
    #Apply the parse on the portion of the datafile containing coordinates
    number,colors,xs,ys,zs = \
        zip(*[parse(line.split()) for line in fraw[iline:iline+natoms]])
    iline+=natoms
    
    cmn=min(colors)
    cmx=max(colors)

    xmn,xmx=min(xs),max(xs)
    ymn,ymx=min(ys),max(ys)
    zmn,zmx=min(zs),max(zs)

    pl.subplot(221)
    pl.scatter(xs,ys,c=colors,vmin=cmn,vmax=cmx,marker='+',faceted=False)
    pl.xlim(xmn,xmx)
    pl.ylim(ymn,ymx)
    pl.title("Z")

    pl.subplot(222)
    pl.scatter(xs,zs,c=colors,vmin=cmn,vmax=cmx,marker='+',faceted=False)
    pl.title("Y")
    pl.xlim(xmn,xmx)
    pl.ylim(zmn,zmx)

    pl.subplot(223)
    pl.scatter(ys,zs,c=colors,vmin=cmn,vmax=cmx,marker='+',faceted=False)
    pl.title("X")
    pl.xlim(ymn,ymx)
    pl.ylim(zmn,zmx)
    
    pl.subplot(224)
    pl.text(0,0.1,info)
    pl.gca().set_xticks([])
    pl.gca().set_yticks([])
    pl.axis('off')

    #pr.prshow("atomProjection.png")
    pl.show()
    pl.gca().clear()
