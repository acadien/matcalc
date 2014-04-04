#!/usr/bin/python

#The ELFCARs analyzed by this script need to be part of rectangular simulations
#parallepiped simulations can be changed to cubic/rectangular using the script
#poscarRectify.py
#This greatly simplifies the application of periodic boundary conditions

import sys
import operator
from scipy import *
from scipy import ndimage
#mine
from struct_tools import minImageDist,minImageAtom,rotmatx
import elfcarIO
from neighbors import voronoiNeighbors
import neighborIO

def elfcarNeighborAnalysis(elfcarfile,verbose=False,minELF=0.5,maxBondLength=4.0):
    #Loads up an ELFCAR, generates a starting neighbor list from voronoi tesselation
    #Initial neighbors are dropped if bondlength is greater than maxBondLength
    #Generates a bunch of points on a cylinder
    #Maps the cylinder points across each neighbor pair
    #If the average ELF is always above minELF accross this cylinder, atoms are bonded.

    #Parse ELFCAR
    elfcar=open(elfcarfile,"r").readlines()
    (basis,atypes,atoms,header),elf = elfcarIO.read(elfcar)
    lengths=array([basis[0][0],basis[1][1],basis[2][2]])

    #Neighbors
    bounds=[[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]
    halfNeighbs = voronoiNeighbors(atoms,basis,style="half")

    a=elf.shape
    AvgElf=sum([sum([sum(line) for line in plane]) for plane in elf])/a[0]/a[1]/a[2]
    if verbose:
        print elfcarfile
        print "Average ELF value:",AvgElf

    #Evaluate the ELF between each nieghbor pair
    #creates a bunch of points on a circle in the x-y plane
    def circlePoints(center,radius):
        N=4
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

    #Loop over atom pairs, generate cylinders
    cylSlices = 10   
    cylRadius = 0.7
    cirPoints = circlePoints(a,cylRadius).T
    cylPoints = [vstack([cirPoints,zeros(cirPoints.shape[1])+float(z)]) for z in range(cylSlices+1)]
    coordination=zeros(len(atoms))
    neighborsELF=[list() for i in range(len(atoms))]

    for i,ineighbs in enumerate(halfNeighbs):
        a=atoms[i]

        for j in ineighbs:
            #copy the cylinder points so you can play with them...
            localCyl=array(cylPoints)

            #Find the minimum image atom
            b=minImageAtom(a,atoms[j],basis)
            l=minImageDist(a,b,basis)

            if l>maxBondLength or l==0.0: continue

            #How to map your cylinder onto the local atom
            R=rotmatx([0,0,l],b-a)

            #Loop over each cylinder slice
            sliceParam=list()
            for cir in localCyl:
                #cir is stretched and rotated cylinder
                cir[2]/=cylSlices/l
                cir=asarray([dot(R,p)+a for p in cir.T]).T
                
                #cir2 is cir with periodic boundary conditions applied and mapped to index space
                cir2=asarray([((c/lengths[ind]*elf.shape[ind])%(elf.shape[ind]-1)) for ind,c in enumerate(cir)])            

                #Interpolation across cylinder
                z=ndimage.map_coordinates(elf,cir2)

                #Average ELF value across each slice.
                sliceParam.append(sum(z)/len(z))

            if min(sliceParam)>minELF:
                #neighbor pair is bonded, count it!
                neighborsELF[i].append(j)
                neighborsELF[j].append(i)

    return neighborsELF

if __name__=="__main__":
    def usage():
        print "Usage:"
        print "%s <ELFCARfile ELFCARfile...> "%(sys.argv[0])
        print "Note, ELFCAR must be part of a cubic/rectangular simulation"

    if len(sys.argv)<2:
        usage()
        exit(0)

    elfcarfiles=sys.argv[1:]

    for elfcarfile in elfcarfiles:
        minELF=0.5
        neighbors = \
            elfcarNeighborAnalysis(elfcarfile,minELF=minELF,maxBondLength=3.8,verbose=True)

        #Write the coordination numbers
        coords=map(len,neighbors)
        data=["Coordination\n"]
        data+=map(lambda x:"%d\n"%(x),coords)
        open("%s.Coord"%elfcarfile,"w").writelines(data)
        
        #Write the coordination histogram
        bins=range(min(coords),max(coords)+1)
        values=[coords.count(i) for i in bins]
        print bins
        print values
        data=["Coordination Count\n"]
        data+=map(lambda x:"%d %d\n"%(x[0],x[1]),zip(bins,values))
        open("%s.CoordHist"%elfcarfile,"w").writelines(data)

        #Write the neighbor list
        #data = ["%d\n"%len(neighbors)]
        #data+= [" ".join(map(str,neighb))+"\n" for neighb in neighbors]
        neighborIO.write("%s.Neighb"%elfcarfile,neighbors)
        #All done.
        print "%s.CoordHist %s.Coord %s.Neighb\n"%(elfcarfile,elfcarfile,elfcarfile)
