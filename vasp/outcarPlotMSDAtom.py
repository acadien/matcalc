#!/usr/bin/python

import sys
from math import *
from scipy import array,zeros
import pylab as pl
#mine
from meanSquareDist import meanSquareDistRefAtom
import outcarIO
#msdIO

#Calculates mean squared distance

def usage():
    print "Usage: %s <Outcar> <Reference Config #>"%sys.argv[0].split("/")[-1]

def outcarMeanSquareDisplaceAtom(outcarFile,refStructure=None):
    outcar=open(outcarFile,"r")
    atoms=list() #atoms[time][atom]

    #Grab ion types
    while True:
        line=outcar.readline()
        if "ions per type" in line:
            break
    atypes=map(int,line.split("=")[1].split())
    Natom=sum(atypes)

    #Grab basis vectors
    while True:
        line=outcar.readline()
        if "direct lattice vectors" in line:
            break
    basis=array([map(float,outcar.readline().split()[:3]) for i in range(3)])
    lengths=array([basis[0][0],basis[1][1],basis[2][2]])

    #Grab atom positions
    count=0
    posit=False
    for line in outcar:
        if posit:
            if "--" in line:
                if len(a)==0:
                    continue
                else:
                    #Analysis
                    atoms.append(array(a))
                    posit=False
            else:
                a.append(map(float,line.split()[:3]))
        elif "POSITION" in line:
            a=list()
            posit=True
            count+=1
    atoms=array(atoms)
    Ntime=len(atoms)

    if refStructure==None:
        delT,msdPerAtom=meanSquareDistRefAtom(atoms,0,Natom,Ntime,lengths)
    elif refStructure > Ntime:
        print "%d requested but %d structures max"%(refStructure,Ntime)
    else:
        delT,msdPerAtom=meanSquareDistRefAtom(atoms,refStructure,Natom,Ntime,lengths)

    return delT,msdPerAtom


if __name__=="__main__":
    if len(sys.argv)<2:
        usage()
        exit(0)
    
    outcarFile=sys.argv[1]
    try:
        num=sys.argv[1].split("_")[1]
    except IndexError:
        num=None

    refStructure=None
    if len(sys.argv)==3:
        refStructure=int(sys.argv[2])

    delT,msdAtom=outcarMeanSquareDisplaceAtom(outcarFile,refStructure)

    pl.figure()
    for atom in range(msdAtom.shape[0]):
        pl.plot(delT,msdAtom[atom])


    dummy,dummy,basis,atoms,forces,types=outcarIO.outcarReadConfig(outcarFile,-1)
    from mayavi import mlab
    
    ax,ay,az=atoms.T
    ops=msdAtom.T[-1]
    v1,v2,v3=basis

    #Box... bleh I hate this code.
    mlab.plot3d([0,v1[0],v1[0]+v2[0],v2[0],0,v3[0]],[0,v1[1],v1[1]+v2[1],v2[1],0,v3[1]],[0,v1[2],v1[2]+v2[2],v2[2],0,v3[2]],color=(0,0,0),line_width=0.5)
    mlab.plot3d([v3[0],v3[0]+v1[0],v3[0]+v2[0]+v1[0],v3[0]+v2[0],v3[0]],[v3[1],v3[1]+v1[1],v3[1]+v2[1]+v1[1],v3[1]+v2[1],v3[1]],[v3[2],v3[2]+v1[2],v3[2]+v2[2]+v1[2],v3[2]+v2[2],v3[2]],color=(0,0,0),line_width=0.5)
    mlab.plot3d([v1[0],v1[0]+v3[0]],[v1[1],v1[1]+v3[1]],[v1[2],v1[2]+v3[2]],color=(0,0,0),line_width=0.5)
    mlab.plot3d([v2[0],v2[0]+v3[0]],[v2[1],v2[1]+v3[1]],[v2[2],v2[2]+v3[2]],color=(0,0,0),line_width=0.5)
    mlab.plot3d([v1[0]+v2[0],v1[0]+v2[0]+v3[0]],[v1[1]+v2[1],v1[1]+v2[1]+v3[1]],[v1[2]+v2[2],v1[2]+v2[2]+v3[2]],color=(0,0,0),line_width=0.5)

    #Atoms
    mp3d = mlab.points3d(ax,ay,az,ops,colormap='gist_earth',scale_factor=1.9,scale_mode='none',resolution=25)
    cb = mlab.colorbar(title=r"MSD ($\AA^2$)", orientation='vertical')#, nb_labels=n,nb_colors=n)

    mn=0
    mx=ceil(max(ops))
    cb.data_range = (mn,mx)
    
    #To make a pretty "hi res" figure use...
    #f=mlab.gcf()
    #f.scene.render_window.aa_frames = 8
    #mlab.draw() 
    pl.show()
