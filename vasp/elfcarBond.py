#!/usr/bin/python

#The ELFCARs analyzed by this script need to be part of rectangular simulations
#parallepiped simulations can be changed to cubic/rectangular using the script
#poscarRectify.py
#This greatly simplifies the application of periodic boundary conditions

import sys
import operator
from scipy import *
from scipy import ndimage,interpolate
import pylab as pl
import mpl_toolkits.mplot3d.axes3d as p3
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import pyplot as plt
#mine
from struct_tools import minImageDist,minImageAtom,rotmatx
import elfcarIO
from neighbors import voronoiNeighbors,full2half
import colors

def elfcarBondAnalysis(elfcarFile,elfcarNeighbsFile,verbose=True):
    #Loads up an ELFCAR, generates a starting neighbor list from voronoi tesselation
    #Initial neighbors are dropped if bondlength is greater than maxBondLength
    #Generates a bunch of points on a cylinder
    #Maps the cylinder points across each neighbor pair
    #If the average ELF is always above minELF accross this cylinder, atoms are bonded.

    ELFlevel=0.5

    #Parse ELFCAR
    elfcar=open(elfcarFile,"r").readlines()
    (basis,atypes,atoms,header),elf = elfcarIO.read(elfcar)
    lengths=array([basis[0][0],basis[1][1],basis[2][2]])

    #Neighbors
    elfcarNeighbs=open(elfcarNeighbsFile,"r").readlines()
    neighbors=[map(int,line.split()) for line in elfcarNeighbs[1:]]
    halfNeighbs=full2half(neighbors)

    a=elf.shape
    AvgElf=sum([sum([sum(line) for line in plane]) for plane in elf])/a[0]/a[1]/a[2]
    if verbose:
        print elfcarFile
        print "Average ELF value:",AvgElf

    #Evaluate the ELF between each nieghbor pair
    #creates a bunch of points on a circle in the x-y plane
    def circlePoints(center,radius):
        N=11
        ps=[float(p)/N*radius for p in range(-N,N+1)]
        r2=radius*radius
        points=list()
        for x in ps:
            for y in ps:
                if (x**2+y**2)<=r2:
                    points.append(asarray([x,y]))
        return asarray(points)

    atoms=asarray([lengths[i]-v for i,v in enumerate(atoms.T)]).T

    #Flip the ELF because its off, thanks a lot VASP.
    elf=elf[::-1,::-1,::-1]

    #Loop over atom pairs
    cylRadius = 2.0
    cirPoints = circlePoints(a,cylRadius).T
    cirPoints = vstack([cirPoints,zeros(cirPoints.shape[1])]) 
    coordination=zeros(len(atoms))
    gridx,gridy=np.mgrid[0:1:20j,0:1:20j]*4.2-2.1
    cirPoints2=array(cirPoints)
    cirPoints2=cirPoints2[:2].T
    f,ax  =plt.subplots(1,1)
    rs=[list() for i in halfNeighbs]
    rsFlat=list()
    for i,ineighbs in enumerate(halfNeighbs):
        a=atoms[i]
        for j in ineighbs:
            #copy the cylinder points so you can play with them...
            localCir=array(cirPoints)

            #Find the minimum image atom
            b=minImageAtom(a,atoms[j],basis)
            l=minImageDist(a,b,basis)

            #How to map your cylinder onto the local atom
            R=rotmatx([0,0,l],b-a)

            #Loop over each cylinder slice
            sliceParam=list()
            localCir[2]=l/2.
            cir=asarray([dot(R,p)+a for p in localCir.T]).T
                
            #cir2 is circle with periodic boundary conditions applied and mapped to index space
            cir2=asarray([((c/lengths[ind]*elf.shape[ind])%(elf.shape[ind]-1)) for ind,c in enumerate(cir)])            

            #Interpolation across cylinder
            z=ndimage.map_coordinates(elf,cir2)

            #Generate a grid and find the desired contour on the bond crosssection
            gridz=interpolate.griddata(cirPoints2,z,(gridx,gridy),fill_value=0)
            ac=ax.contour(gridx,gridy,gridz,levels=[ELFlevel],linewidths=5,colors="black")

            #take the longest list (this is the bond)
            m=0
            if len(ac.collections[0].get_paths()):
                m=np.asarray([len(v.vertices) for v in ac.collections[0].get_paths()]).argmax()     
            midCircle= ac.collections[0].get_paths()[m].vertices.T
            n=midCircle.shape[1]

            #Find the center of the bond and calculate the radius
            com=[midCircle[0].sum()/n,midCircle[1].sum()/n]
            midCircle[0]-=com[0]
            midCircle[1]-=com[1]
            r=0
            for m in midCircle.T:
                r+=(m[0]**2+m[1]**2)**(0.5)
            r/=n
            rs[i].append(r)
            rs[j].append(r)
            rsFlat.append(r)

#                pl.contour(gridx,gridy,gridz,levels=[0.5],linewidths=5,colors="black")
#                pl.show()
            #avgCircle+=jCircle
    #avgCircle/=sum(map(len,neighbors))/2.
    
    #pl.show()
    return rs,rsFlat,atoms,neighbors,basis

if __name__=="__main__":
    def usage():
        print "Usage:"
        print "%s <ELFCARfile ELFCARfile...> "%(sys.argv[0])
        print "Loads ELFCAR and ELFCAR.neighbors, does analysis of bonds between neighbors."


    if len(sys.argv)<2:
        usage()
        exit(0)

    elfcarfiles=sys.argv[1:]

    for elfcarfile in elfcarfiles:
        ELFlevel=0.5
        #radius of bonds
        rs,rsFlat,atoms,neighbors,basis=elfcarBondAnalysis(elfcarfile,elfcarfile+".Neighb",verbose=True)
    
    """
    pl.show()
    pl.clf()
    pl.hist(rsFlat,bins=[i/35.+0.5 for i in range(100)])
    pl.show()
    """

        #Plotting of atoms with bond radius.
        from mayavi import mlab
        fig=mlab.figure(bgcolor=(1.0,1.0,1.0))

        #coordination of atoms
        coord=map(lambda x: float(x.split()[0]), open(elfcarfile+".Coord","r").readlines()[1:])
        mx=max(coord)
        mn=min(coord)
        dl=mx-mn
        coord=map(lambda x: (x-mn)/dl,coord)
        mlab.points3d(atoms[:,0],atoms[:,1],atoms[:,2],coord,scale_factor=1.0,scale_mode='none')

        #bond radius
        mx=max(rsFlat)
        mn=min(rsFlat)
        dl=mx-mn
        for i,ineighbs in enumerate(neighbors):
            a=atoms[i]
            for i2,j in enumerate(ineighbs):
                b=minImageAtom(a,atoms[j],basis)
                c=(rs[i][i2]-mn)/dl
                mlab.plot3d([a[0],b[0]],[a[1],b[1]],[a[2],b[2]],color=(c,1-c,1-c),tube_radius=0.1,tube_sides=10)
        mlab.show()

    #need to do this plot separately
    pl.plot([a[0],b[0]],[a[1],b[1]],[a[2],b[2]],c=r/1.6,lw=10)
