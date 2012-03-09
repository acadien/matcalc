#!/usr/bin/python

#The ELFCAR's analyzed by this script need to be part of rectangular simulations
#parallepiped simulations can be changed to cubic/rectangular using the script
#poscarRectify.py
#This greatly simplifies the application of periodic boundary conditions

import sys
import operator
from scipy import *
#mine
from poscarIO import readposcar
from datatools import flatten
from struct_tools import *
from voronoiNeighbors import voronoiNeighbors

#Checks if the local-grid is within the bounds of the global-grid
def outOfBounds(lgrids,bounds):
    for lgrid in lgrids:
        for i in range(3):
            if lgrid[i][1]>bounds[i]:
                return True
            elif lgrid[i][0]<0:
                return True
    return False

def whichBounds(lgrids,bounds):
    for a,lgrid in enumerate(lgrids):
        for i in range(3):
            if lgrid[i][1]>bounds[i]:
                return a
            elif lgrid[i][0]<0:
                return a

def bilinear_interpolation(x, y, points):
    points = sorted(points)
    (x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points
    return (q11 * (x2 - x) * (y2 - y) +  q21 * (x - x1) * (y2 - y) +  \
            q12 * (x2 - x) * (y - y1) +  q22 * (x - x1) * (y - y1) )  \
            / ((x2 - x1) * (y2 - y1) + 0.0)

def interp3d(ipnt,points,bnds,data):
    [xp,yp,zp]=points
    np=map(lambda x:max(x,0),[sum([1 for j in points[i] if ipnt[i]>j])-1 for i in range(3)])

    axs=[np[0]]*4
    bxs=[np[0]+1]*4
    ys=[np[1],np[1]]+[np[1]+1,np[1]+1]
    zs=[np[2],np[2]+1]*2

    #Use bilinear interpolation to get values on opposing planes
    adata=[data[x][y][z] for x,y,z in zip(axs,ys,zs)]
    aval=bilinear_interpolation(ipnt[1],ipnt[2],zip([points[1][i] for i in ys],[points[2][i] for i in zs],adata))
    bdata=[data[x][y][z] for x,y,z in zip(bxs,ys,zs)]
    bval=bilinear_interpolation(ipnt[1],ipnt[2],zip([points[1][i] for i in ys],[points[2][i] for i in zs],bdata))

    #Weighted average of planar values
    return (aval*fabs(ipnt[0]-xp[bxs[0]])+bval*fabs(xp[axs[0]]-ipnt[0]))/fabs(xp[bxs[0]]-xp[axs[0]])


def linePoints(a,b,N): #Creates a 3D line from 2 points
    eqn=lambda t: a+t*(b-a)
    return [eqn(float(i)/N) for i in range(N)]

def elfNeighbors(atoms,basis,atypes,gridSize,elf,halfNeighbors=None,bondCutoffs=0.5,Ninterp=50):
    singlet=False
    if type(bondCutoffs)==type(float()):
        singlet=True
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
    if not singlet:
        avgelfs=zeros([len(bondCutoffs),Ninterp])
        avgfails=zeros([len(bondCutoffs),Ninterp])
        naes=zeros(len(bondCutoffs))
        nafs=zeros(len(bondCutoffs))
        eHNs=[[list() for a in range(len(atoms))] for b in range(len(bondCutoffs))]
    else:
        avgelfs=zeros([1,Ninterp])
        avgfails=zeros([1,Ninterp])
        naes=zeros(1)
        nafs=zeros(1)
        eHNs=[[list() for a in range(len(atoms))] for b in range(1)]
    #Local Grid size

    delGrids=simSize/gridSize
    lGridSize=gridSize
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

        lGridCenter = map(int,atomi/simSize*gridSize)
        lGridBounds = [[lGridCenter[i]-lGridSize[i]/2,lGridCenter[i]+lGridSize[i]/2] for i in range(3)]
        lGridPoints = [map(lambda x: x*delGrids[i],range(lGridCenter[i]-lGridSize[i]/2,lGridCenter[i]+lGridSize[i]/2)) for i in range(3)]

        #Stitch together elf in local elf (lelf) to enforce periodicity, use bounds to simplify
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

        lelf=zeros(lGridSize)
        for lbnd,bnd in zip(lbnds,bnds):
            [[lxa,lxb],[lya,lyb],[lza,lzb]]=lbnd
            [[xa,xb],[ya,yb],[za,zb]]=bnd
            lelf[lxa:lxb, lya:lyb, lza:lzb]=elf[xa:xb, ya:yb, za:zb]

        #Loop over Neighboring atoms within the local grid
        cnt=0
        for jNeighb in jNeighbors:
            atomj=latoms[jNeighb]
            for i,aj in enumerate(atomj):
                if aj < lGridPoints[i][0]:
                    atomj[i]+=simSize[i]
                if aj > lGridPoints[i][-1]:
                    atomj[i]-=simSize[i]

            if  atomj[0]>lGridPoints[0][-1] or atomj[1]>lGridPoints[1][-1] or atomj[2]>lGridPoints[2][-1] or \
                atomj[0]<lGridPoints[0][0]  or atomj[1]<lGridPoints[1][0]  or atomj[2]<lGridPoints[2][0]:
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
            #Interpolate on local ELF grid between these two atoms 75% of time spent here
            pnts2interp=linePoints(array([0,0,0]),atomj-atomi,Ninterp) 

            vals=[interp3d(i+atomi,lGridPoints,lGridBounds,lelf) for i in pnts2interp] #<- most time spent here (~5.8seconds)

            if not singlet:
                for j,bondCutoff in enumerate(bondCutoffs):
                    bb=False
                    for i in vals:
                        if i<bondCutoff:
                            cnt+=1
                            bb=True
                            break
                    if bb:
                        nafs[j]+=1
                        avgfails[j]+=array(vals)                   
                    else:
                        naes[j]+=1
                        eHNs[j][iNeighb].append(jNeighb)
                        avgelfs[j]+=array(vals)
            else:
                bb=False
                for i in vals:
                    if i<bondCutoffs:
                        cnt+=1
                        bb=True
                        break
                if bb:
                    nafs[0]+=1
                    avgfails[0]+=array(vals)                   
                else:
                    naes[0]+=1
                    eHNs[0][iNeighb].append(jNeighb)
                    avgelfs[0]+=array(vals)

        #print "worked!",cnt
    for naf,avgfail,nae,avgelf in zip(nafs,avgfails,naes,avgelfs):
        if naf==0:
            avgfail=zeros(len(avgfail))
        else:
            avgfail/=naf
        if nae==0:
            avgelf=zeros(len(avgelf))
        else:
            avgelf/=nae
    if singlet:
        return eHNs[0],avgelfs[0],avgfails[0]
    return eHNs,avgelfs,avgfails


if __name__=="__main__":
    def usage():
        print "Usage:"
        print "%s <ELFCAR> <bond-ELF-cutoff=0.5>"%(sys.argv[0])
        print "Note ELFCAR must be part of a cubic/rectangular simulation"

    if len(sys.argv)<2:
        usage()
        exit(0)

    bondcutoff=[0.5]
    if len(sys.argv)==3:
        bondcutoff=[float(sys.argv[2])]

    #Parse ELFCAR
    elfcar=open(sys.argv[1],"r").readlines()
    v1,v2,v3,atypes,axs,ays,azs,header,elfcar = readposcar(elfcar)
    basis=asarray([v1,v2,v3])
    bounds=[[0.,v1[0]],[0.,v2[1]],[0.,v3[2]]]
    atoms=asarray(zip(axs,ays,azs))

    #Read in Grid properties
    elfcar.pop(0)
    gridSize=[int(i) for i in elfcar.pop(0).split()]
    nGridPoints=reduce(operator.mul,gridSize)

    #Read in ELF data
    elf=asarray([float(i) for i in ("".join(elfcar)).split()])
    elf.shape=gridSize

    halfNeighbors=voronoiNeighbors(atoms=atoms,basis=basis,atypes=atypes,style='half')

    print "Number of Neighbors before elimination:",sum([len(i) for i in halfNeighbors])

    #halfNeighbors,avgelf,avgfail=elfNeighbors(atoms,basis,gridSize,elf,r,bondCutoff=bondcutoff)

    a=elf.shape
    print "Average ELF value:",sum([sum([sum(line)/a[2] for line in plane])/a[1] for plane in elf])/a[0]
        
    fige=pl.figure()
    ae=fige.add_subplot(211)
    af=fige.add_subplot(212)
    eHalfNeighbors,avges,avgfs=elfNeighbors(atoms,basis,atypes,gridSize,elf,halfNeighbors,bondCutoffs=[float(i)/10 for i in range(2,7)]) 
    print "Neighbors after elim:",sum(map(len,eHalfNeighbors))
    for i,avge,avgf in zip(range(2,7),avges,avgfs):
        ae.plot(avge,label=str(float(i)/10.))
        af.plot(avgf)

    #pl.legend(["Good Neighbs","Bad Neighbs"])
    ae.set_title("Average ELF value for Bonds")
    ae.legend(title="cutoff",bbox_to_anchor=(1.10,1,))
    af.set_title("Average ELF value for non-Bonds")

    pl.show()
    exit(0)

    
