#!/usr/bin/python

import sys,operator
from scipy import *
from scipy.interpolate import griddata
#mine
from struct_tools import dist_periodic,dist,rotmatx,mag
from interpolate import interp1d,interp3d

#These functions interpolate 3D gridded fields about a set of points
#Particularily useful for working with CHGCARs and neighbors.

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


#Creates a 3D grid inside the box using the number of points specified
#box=[Xlow,Xhigh,Ylow,Yhigh,Zlow,Zhigh]
#Ns=[Nx,Ny,Nz]
def volumePoints(box,Ns):
    Nx,Ny,Nz=Ns
    X=zeros(Nx)
    Y=zeros(Ny)
    Z=zeros(Nz)
    
    dx=box[1]-box[0]
    X=[box[0]+i*dx for i in range(Nx)]
    dy=box[3]-box[2]
    Y=[box[2]+i*dy for i in range(Ny)]
    dz=box[5]-box[4]
    Z=[box[4]+i*dz for i in range(Nz)]
    
    return X,Y,Z


#Loop over neighbors provided (otherwise use vornoiNeighbors)
#Interpolate the 3D Field between each of these neighbors
def fieldNeighbors1D(atoms,atypes,basis,field,halfNeighbors=None,Ninterp=50):
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
    fieldSize=elf.size
    delGrids=simSize/fieldSize
    lGridSize=fieldSize
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

        lGridCenter = map(int,atomi/simSize*fieldSize)
        lGridBounds = [[lGridCenter[i]-lGridSize[i]/2,lGridCenter[i]+lGridSize[i]/2] for i in range(3)]
        lGridBoundPoints = array([[lGridBounds[i][0]*delGrids[i],lGridBounds[i][1]*delGrids[i]] for i in range(3)])
        lGridPoints = array([map(lambda x: x*delGrids[i],range(lGridCenter[i]-lGridSize[i]/2,lGridCenter[i]+lGridSize[i]/2)) for i in range(3)])

        #Stitch together field in local field (lfield) to enforce periodicity, use bounds to simplify
        bnds=list()
        lbnds=list()
        bnds.append(list(lGridBounds))
        lbnds.append(map(list,zip([0,0,0],lGridSize)))
        while outOfBounds(bnds,fieldSize):
            a=whichBounds(bnds,fieldSize)
            bnds.append(map(list,bnds[a]))
            lbnds.append(map(list,lbnds[a]))
            for i in range(3):
                if bnds[-1][i][1] > fieldSize[i]:
                    lbnds[a][i][1] = lGridSize[i]-(bnds[-1][i][1]-fieldSize[i])
                    lbnds[-1][i][0] = lGridSize[i]-(bnds[-1][i][1]-fieldSize[i])
                    lbnds[-1][i][1] = lbnds[-1][i][0]+(bnds[-1][i][1]-fieldSize[i])

                    bnds[-1][i][0]  = 0
                    bnds[-1][i][1] -= fieldSize[i]
                    bnds[a][i][1]   = fieldSize[i]
                    break

                elif bnds[-1][i][0] < 0:
                    lbnds[-1][i][0] = 0
                    lbnds[-1][i][1] = -bnds[-1][i][0]
                    lbnds[a][i][0] = -bnds[-1][i][0]

                    bnds[-1][i][1]  = fieldSize[i]
                    bnds[-1][i][0] += fieldSize[i]
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
            pnts2interp=linePoints(atomi,atomj,Ninterp) 
            yys=array([interp3d(i,lGridBoundPoints,lGridPoints,lfield) for i in pnts2interp])
            xwid=dist(pnts2interp[0],pnts2interp[-1])
            xdel=xwid/Ninterp
            xxs=[xdel*x-xwid/2 for x in range(Ninterp)]
            xlines[-1][-1]=xxs
            ylines[-1][-1]=yys
            avgyline+=yys

    avgyline/=nlines
    return avgyline,xlines,ylines,halfNeighbors

#Loop over neighbors provided (otherwise use vornoiNeighbors)
#Interpolate the 3D Field between each of these neighbors on to a 3D Grid about the atom pairs
#Note: it is much more efficient to do the bondcutoff comparision within the grid loop, as weird as it may be to work with.
def fieldNeighbors3D(atoms,atypes,basis,field,fieldSize,halfNeighbors=None,Ninterps=[7,7,15],cutoffs=None,loc=None):
    #loc=["center","bond","half"]

    [v1,v2,v3]=basis
    latoms=array(atoms)

    if loc==None:
        loc="bond"

    #Ensure simulation box is orthorhombic
    if sum(map(fabs,v1))-fabs(v1[0]) > 1e-10 or \
       sum(map(fabs,v2))-fabs(v2[1]) > 1e-10 or \
       sum(map(fabs,v3))-fabs(v3[2]) > 1e-10:
        print "Error: simulation axes are not sufficiently orthogonal.  Use script poscarRectify.py and regenerate doscar"
        exit(0)
    
    cutFlag=True
    if cutoffs==None:
        cutFlag=False
    else: #ensure cutoffs is a list and not just a float
        cutoffs=list(cutoffs)

    bounds=[[0.,v1[0]],[0.,v2[1]],[0.,v3[2]]]
    simSize=asarray([v1[0],v2[1],v3[2]])

    if halfNeighbors==None:
        halfNeighbors=voronoiNeighbors(atoms=latoms,basis=basis,atypes=atypes,style='half')

    #Interpolation and summation properties
    grids=list()
    if cutFlag:
        avgGrids=[zeros(Ninterps) for i in range(len(cutoffs))]
        bondcnts=zeros(len(cutoffs))
    else:
        avgGrid=zeros(Ninterps)
        bondcnt=0

    nlines=float(sum(map(len,halfNeighbors)))
    delx=2.0/Ninterps[0]
    dely=2.0/Ninterps[1]

    #Local Grid size
    delGrids=simSize/fieldSize
    lGridSize=fieldSize
    nlGridPoints=reduce(operator.mul,lGridSize)
    for iNeighb,jNeighbors in enumerate(halfNeighbors):
        if len(jNeighbors)==0:
            continue

        atomi=latoms[iNeighb]
        for i,ai in enumerate(atomi):
            if ai < 0.0:
                atomi[i]+=simSize[i]
            if ai > simSize[i]:
                atomi[i]-=simSize[i]
        
        #Make the local grid about the central atom
        lGridCenter = map(int,atomi/simSize*fieldSize)
        lGridBounds = [[lGridCenter[i]-lGridSize[i]/2,lGridCenter[i]+lGridSize[i]/2] for i in range(3)]
        lGridBoundPoints = array([[lGridBounds[i][0]*delGrids[i],lGridBounds[i][1]*delGrids[i]] for i in range(3)])
        lGridPoints = array( [ \
                map(lambda x: x*delGrids[i],range(lGridCenter[i]-lGridSize[i]/2,lGridCenter[i]+lGridSize[i]/2)) \
                    for i in range(3) ] )
        bnds=list()
        lbnds=list()
        bnds.append(list(lGridBounds))
        lbnds.append(map(list,zip([0,0,0],lGridSize)))
        while outOfBounds(bnds,fieldSize):
            a=whichBounds(bnds,fieldSize)
            bnds.append(map(list,bnds[a]))
            lbnds.append(map(list,lbnds[a]))
            for i in range(3):
                if bnds[-1][i][1] > fieldSize[i]:
                    lbnds[a][i][1] = lGridSize[i]-(bnds[-1][i][1]-fieldSize[i])
                    lbnds[-1][i][0] = lGridSize[i]-(bnds[-1][i][1]-fieldSize[i])
                    lbnds[-1][i][1] = lbnds[-1][i][0]+(bnds[-1][i][1]-fieldSize[i])

                    bnds[-1][i][0]  = 0
                    bnds[-1][i][1] -= fieldSize[i]
                    bnds[a][i][1]   = fieldSize[i]
                    break

                elif bnds[-1][i][0] < 0:
                    lbnds[-1][i][0] = 0
                    lbnds[-1][i][1] = -bnds[-1][i][0]
                    lbnds[a][i][0] = -bnds[-1][i][0]

                    bnds[-1][i][1]  = fieldSize[i]
                    bnds[-1][i][0] += fieldSize[i]
                    bnds[a][i][0]   = 0
                    break

        #Now make the field
        lfield=zeros(lGridSize)
        for lbnd,bnd in zip(lbnds,bnds):
            [[lxa,lxb],[lya,lyb],[lza,lzb]]=lbnd
            [[xa,xb],[ya,yb],[za,zb]]=bnd
            lfield[lxa:lxb, lya:lyb, lza:lzb]=field[xa:xb, ya:yb, za:zb]

        #Loop over Neighboring atoms within the local grid
        for jNeighb in jNeighbors:
            atomj=latoms[jNeighb]
            for i,aj in enumerate(atomj):
                if aj < lGridPoints[i][0]:
                    atomj[i]+=simSize[i]
                if aj > lGridPoints[i][-1]:
                    atomj[i]-=simSize[i]
            d=dist(atomj,atomi)
            
            #If this bond is longer than all the cutoffs than don't bother with it
            if cutFlag and len([1 for i in cutoffs if d<=i])==0:
                continue

            #Change the location as requested
            if loc=="bond":
                delz=mag(atomj-atomi)/Ninterps[2]
                Tran=(atomj-atomi)/2.+atomi
            if loc=="half":
                delz=mag(atomj-atomi)/(Ninterps[2]*2.)
                Tran=(atomj-atomi)/4.+atomi
            if loc=="center":
                delz=dely
                Tran=atomi

            #Change of coordinate systems into atomj-atomi basis.
            Rota=rotmatx(array([0.,0.,1.]),atomj-atomi) 

            xps=array([(i-Ninterps[0]/2)*delx for i in range(Ninterps[0])]*Ninterps[1]*Ninterps[2]).ravel()
            yps=array([[(i-Ninterps[1]/2)*dely]*Ninterps[0] for i in range(Ninterps[1])]*Ninterps[2]).ravel()
            zps=array([[(i-Ninterps[2]/2)*delz]*Ninterps[0]*Ninterps[1] for i in range(Ninterps[2])]).ravel()
            ipnts=array([dot(Rota,array(pnt))+Tran for pnt in zip(xps,yps,zps)])

            #Do the interpolation
            if cutFlag:
                for j,cut in enumerate(cutoffs):
                    if d <= cut: 
                        bondcnts[j]+=1
                        if Ninterps[0]!=1:
                            grid=array([interp3d(ipnts[i],lGridBoundPoints,lGridPoints,lfield) for i in range(len(ipnts))]).reshape(Ninterps)
                        else:
                            grid=array([interp3d(ipnts[i],lGridBoundPoints,lGridPoints,lfield) for i in range(len(ipnts))]).reshape(Ninterps[1:])
                        grids.append(grid)
                        avgGrids[j]+=grid
            else:    
                grid=array([interp3d(ipnts[i],lGridBoundPoints,lGridPoints,lfield) for i in range(len(ipnts))]).reshape(Ninterps)
                grids.append(grid)
                avgGrid+=grid
                bondcnt+=1
    if cutFlag:
        for i,v in enumerate(bondcnts):
            if v==0:
                bondcnts[i]=1
        return [avgGrid/bc for bc,avgGrid in zip(bondcnts,avgGrids)],grids,bondcnts
    return avgGrid/bondcnt,grids,bondcnt
