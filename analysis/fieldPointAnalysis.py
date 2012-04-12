#!/usr/bin/python

import sys,operator
from scipy import *
from scipy.interpolate import griddata
#mine
from struct_tools import dist_periodic,dist,rotmatx,mag
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
#Note: it is much more efficient to do the bondcutoff comparision within the grid loop, as weird as it may be to work with.
def fieldNeighbors1D(atoms,atypes,basis,field,fieldSize,halfNeighbors=None,Ninterp=50):
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
#Interpolate the 3D Field between each of these neighbors on to a 3D Grid
#Note: it is much more efficient to do the bondcutoff comparision within the grid loop, as weird as it may be to work with.
def fieldNeighbors3D(atoms,atypes,basis,field,fieldSize,halfNeighbors=None,Ninterps=[5,5,20]):
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
    avgGrid=zeros(Ninterps)
    nlines=float(sum(map(len,halfNeighbors)))
    delx=1.0/Ninterps[0]
    dely=1.0/Ninterps[1]


    #Local Grid size
    delGrids=simSize/fieldSize
    lGridSize=fieldSize
    nlGridPoints=reduce(operator.mul,lGridSize)
    for iNeighb,jNeighbors in enumerate(halfNeighbors):
        #xlines.append(list())
        #ylines.append(list())
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

#        xps=array([i for i in lGridPoints[0]]*lGridSize[1]*lGridSize[2]).ravel()
#        yps=array([[i]*lGridSize[0] for i in lGridPoints[1]]*lGridSize[2]).ravel()
#        zps=array([[i]*lGridSize[0]*lGridSize[1] for i in lGridPoints[2]]).ravel()
#        lpoints=array(zip(xps,yps,zps))
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
            #xlines[-1].append(zeros([Ninterp]))
            #ylines[-1].append(zeros([Ninterp]))
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
            #Interpolate on local FIELD grid between these two atoms... 75% of time spent here
            #pnts2interp=linePoints(array([0,0,0]),atomj-atomi,Ninterp) 
            R=rotmatx(array([0,0,1]),atomj-atomi) 
            T=(atomj-atomi)/2+atomi
            delz=mag(atomj-atomi)/Ninterps[2]

            #Change of coordinate systems into atomj-atomi basis.
            #Translate so center of cube is at 0.
            #Perform rotation on points
            #Translate back.
            xps=array([(i-Ninterps[0]/2+0.5)*delx for i in range(Ninterps[0])]*Ninterps[1]*Ninterps[2]).ravel()
            yps=array([[(i-Ninterps[1]/2+0.5)*dely]*Ninterps[0] for i in range(Ninterps[1])]*Ninterps[2]).ravel()
            zps=array([[(i-Ninterps[2]/2+0.5)*delz]*Ninterps[0]*Ninterps[1] for i in range(Ninterps[2])]).ravel()
            ipnts=array([dot(R,array(pnt))+T for pnt in zip(xps,yps,zps)])
            #xps,yps,zps=zip(*ipnts)
            #lGridPoints=lGridPoints.T#array(zip(*lGridPoints))
            #lGridPoints.shape=len(lGridPoints),3
            #print lGridPoints.shape
            #print ipnts.shape
            avgGrid+=array([interp3d(ipnts[i],lGridBoundPoints,lGridPoints,lfield) for i in range(len(ipnts))]).reshape(Ninterps)
            #print avgGrid    
            #avgGrid+=griddata(lGridPoints,lfield,ipnts)
            #print avgGrid
            
            """
            print min(xps),max(xps)
            print min(yps),max(yps)
            print min(zps),max(zps)
            import matplotlib as mpl
            from mpl_toolkits.mplot3d import Axes3D
            import matplotlib.pyplot as plt

            fig = plt.figure()
            ax = fig.gca(projection='3d')
            print len(lpoints)
            for i in range(0,len(lpoints),1500):
                ax.scatter(*(lpoints[i]),c="blue")
            for i in range(Ninterps[0]*Ninterps[1]*Ninterps[2]):
                ax.scatter(xps[i],yps[i],zps[i],c="red")
            plt.show()
            #ax.scatter(*list(xps[1]),c="red")
            #ax.scatter(*list(xps[2]),c="red")
            #ax.scatter(*list(xps[3]),c="red")

            #ax.scatter(*atomj)
            #ax.scatter(*atomi)

            #plt.show()
            #spec=lambda x:"array(["+",".join(map(str,x))+"])"
            #print "xs="+spec(xps)
            #print "ys="+spec([spec(i) for i in yps])
            #print "zs="+spec([spec(i) for i in zps])
                


            
            cx=atomj[x
            wx=wy=3.0
            boxi=[,x,0,y,0,z]
            Xi,Yi,Zi=volumePoints(
                yys=array([interp3d(i+atomi,lGridBoundPoints,lGridPoints,lfield) for i in pnts2interp])
            xwid=dist(pnts2interp[0],pnts2interp[-1])
            xdel=xwid/Ninterp
            xxs=[xdel*x-xwid/2 for x in range(Ninterp)]
            #xlines[-1][-1]=xxs
            #ylines[-1][-1]=yys
            avgyline+=yys
            """
    avgyline/=nlines
    return avgGrid,halfNeighbors
