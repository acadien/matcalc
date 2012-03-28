#!/usr/bin/python

import sys,operator
from scipy import *
#mine
from struct_tools import dist_periodic,dist
from interpolate import interp1d,interp3d

#delete me
import pylab as pl

#Returns True if any point in lgrids is out of 'bounds'
def outOfBounds(lgrids,bounds):
    for lgrid in lgrids:
        for i in range(3):
            if lgrid[i][1]>bounds[i]:
                return True
            elif lgrid[i][0]<0:
                return True
    return False

#Returns the index of the point that is out of bounds
def whichBounds(lgrids,bounds):
    for a,lgrid in enumerate(lgrids):
        for i in range(3):
            if lgrid[i][1]>bounds[i]:
                return a
            elif lgrid[i][0]<0:
                return a

#Creates a 3D line from 2 points
def linePoints(a,b,N): 
    eqn=lambda t: a+t*(b-a)
    return [eqn(float(i)/N) for i in range(N)]

def addline(xx,yy,xadd,yadd):
    yyadd=interp1d(xadd,yadd,xx)
    return yy+yyadd

#Loop over neighbors provided (otherwise use vornoiNeighbors)
#Interpolate the 3D Field between each of these neighbors
#Note: it is much more efficient to do the bondcutoff comparision within the grid loop, as weird as it may be to work with.
def fieldNeighbors(atoms,atypes,basis,field,gridSize,halfNeighbors=None,Ninterp=50):
    [v1,v2,v3]=basis
    latoms=array(atoms)

    #Ensure simulation box is orthorhombic
    if sum(map(fabs,v1))-fabs(v1[0]) > 1e-10 or \
       sum(map(fabs,v2))-fabs(v2[1]) > 1e-10 or \
       sum(map(fabs,v3))-fabs(v3[2]) > 1e-10:
        print "Error: simulation axes are not sufficiently orthogonal.  Use script poscarRectify.py and regenerate doscar"
        exit(0)
    
    bounds=[[0.,v1[0]],[0.,v2[1]],[0.,v3[2]]]
    simSize=asarray([v1[0],v2[1],v3[2]])

    if halfNeighbors==None:
        halfNeighbors=voronoiNeighbors(atoms=latoms,basis=basis,atypes=atypes,style='half')

    #Interpolation and summation properties
    xlines=list()
    ylines=list()
    avgyline=zeros(Ninterp)
    nlines=float(sum(map(len,halfNeighbors)))

    #Local Grid size
    delGrids=simSize/gridSize
    lGridSize=gridSize
    nlGridPoints=reduce(operator.mul,lGridSize)
    for iNeighb,jNeighbors in enumerate(halfNeighbors):
        xlines.append(list())
        ylines.append(list())
        if len(jNeighbors)==0:
            continue

        atomi=latoms[iNeighb]
        for i,ai in enumerate(atomi):
            if ai < 0.0:
                atomi[i]+=simSize[i]
            if ai > simSize[i]:
                atomi[i]-=simSize[i]

        lGridCenter = map(int,atomi/simSize*gridSize)
        lGridBounds = [[lGridCenter[i]-lGridSize[i]/2,lGridCenter[i]+lGridSize[i]/2] for i in range(3)]
        lGridBoundPoints = array([[lGridBounds[i][0]*delGrids[i],lGridBounds[i][1]*delGrids[i]] for i in range(3)])
        lGridPoints = array([map(lambda x: x*delGrids[i],range(lGridCenter[i]-lGridSize[i]/2,lGridCenter[i]+lGridSize[i]/2)) for i in range(3)])

        #Stitch together field in local field (lfield) to enforce periodicity, use bounds to simplify
        bnds=list()
        lbnds=list()
        bnds.append(list(lGridBounds))
        lbnds.append(map(list,zip([0,0,0],lGridSize)))
        while outOfBounds(bnds,gridSize):
            a=whichBounds(bnds,gridSize)
            bnds.append(map(list,bnds[a]))
            lbnds.append(map(list,lbnds[a]))
            for i in range(3):
                if bnds[-1][i][1] > gridSize[i]:
                    lbnds[a][i][1] = lGridSize[i]-(bnds[-1][i][1]-gridSize[i])
                    lbnds[-1][i][0] = lGridSize[i]-(bnds[-1][i][1]-gridSize[i])
                    lbnds[-1][i][1] = lbnds[-1][i][0]+(bnds[-1][i][1]-gridSize[i])

                    bnds[-1][i][0]  = 0
                    bnds[-1][i][1] -= gridSize[i]
                    bnds[a][i][1]   = gridSize[i]
                    break

                elif bnds[-1][i][0] < 0:
                    lbnds[-1][i][0] = 0
                    lbnds[-1][i][1] = -bnds[-1][i][0]
                    lbnds[a][i][0] = -bnds[-1][i][0]

                    bnds[-1][i][1]  = gridSize[i]
                    bnds[-1][i][0] += gridSize[i]
                    bnds[a][i][0]   = 0
                    break

        lfield=zeros(lGridSize)
        for lbnd,bnd in zip(lbnds,bnds):
            [[lxa,lxb],[lya,lyb],[lza,lzb]]=lbnd
            [[xa,xb],[ya,yb],[za,zb]]=bnd
            lfield[lxa:lxb, lya:lyb, lza:lzb]=field[xa:xb, ya:yb, za:zb]

        #Loop over Neighboring atoms within the local grid
        for jNeighb in jNeighbors:
            xlines[-1].append(zeros([Ninterp]))
            ylines[-1].append(zeros([Ninterp]))
            atomj=latoms[jNeighb]
            for i,aj in enumerate(atomj):
                if aj < lGridPoints[i][0]:
                    atomj[i]+=simSize[i]
                if aj > lGridPoints[i][-1]:
                    atomj[i]-=simSize[i]

            #Error detection
            if  atomj[0]>lGridPoints[0][-1] or atomj[1]>lGridPoints[1][-1] or atomj[2]>lGridPoints[2][-1] or \
                atomj[0]<lGridPoints[0][0]  or atomj[1]<lGridPoints[1][0]  or atomj[2]<lGridPoints[2][0]:
                d=dist_periodic(atomi,atomj,array(zip([0,0,0],simSize)))
                if d>7.5: continue
                print "="*50
                print "Start Error Report:"
                print "AtomI(%d):"%iNeighb,atomi
                print "AtomJ(%d):"%jNeighb,atomj
                print "Distance:",dist_periodic(atomi,atomj,array(zip([0,0,0],simSize)))
                print "GridX Min,Max:","% 5.5f % 5.5f"%(min(lGridPoints[0]),max(lGridPoints[0]))
                print "GridY Min,Max:","% 5.5f % 5.5f"%(min(lGridPoints[1]),max(lGridPoints[1]))
                print "GridZ Min,Max:","% 5.5f % 5.5f"%(min(lGridPoints[2]),max(lGridPoints[2]))

                print "Error: AtomJ seems to be lying outside the local grid, which should not be possible."
                exit(0)
            
            #This is where the magic happens, finally do some computing...
            #Interpolate on local FIELD grid between these two atoms 75% of time spent here
            pnts2interp=linePoints(array([0,0,0]),atomj-atomi,Ninterp) 
            yys=array([interp3d(i+atomi,lGridBoundPoints,lGridPoints,lfield) for i in pnts2interp])
            xxs=map(lambda x: dist(pnts2interp[0],x),pnts2interp)
            xlines[-1][-1]=xxs
            ylines[-1][-1]=yys
            avgyline+=yys

    avgyline/=nlines
    return avgyline,xlines,ylines,halfNeighbors
