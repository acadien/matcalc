#!/usr/bin/python
import scipy
from numpy import *
#mine
import gaussFunctions

#Helper function finds the header for next Kpoint and returns it at the ion label
def grabKP(procar):
    while True:
        line=procar.readline()
        if "k-point" in line:
            break
    line=line.split("=")
    weight=float(line[1].strip())
    kpoint=line[0].split(":")[1].split()[:-1]
    if len(kpoint) !=3:
        temp=list()
        for i,j in enumerate(map(lambda x: x[1:].find("-"),kpoint)):
            if j==-1:
                temp.append(kpoint[i])
            else:
                temp+=[kpoint[i][:j+1],kpoint[i][j+1:]]
        kpoint=temp
    kpoint=map(float,kpoint)            
    procar.readline()
    return kpoint,weight

def grabE(procar):
    while True:
        line=procar.readline()
        if "energy" in line:
            break
    energy=float(line.split()[4])
    procar.readline()
    return energy

#Helper function grabs the occupancies from the PROCAR and averages over the number of ions for a band
def grabOatB(procar, Nions):
    header=procar.readline().split()[1:]
    Norbs=len(header)

    for i in range(Nions):
        procar.readline()

    #occ=map(lambda x:sum(x),zip(*[map(float,procar.readline().split()[1:]) for i in range(Nions)]))
    occ = array(map(float,procar.readline().split()[1:]))
    return occ,header

#Generate a bunch of gaussians, centered at x with height y and sum the up on the grid sg
#xgrid and ygrid must be arrays
def sumgauss(xcenter,ycenter,xgrid,ygrid,sigma):
    return gaussFunctions.gauss1D(xgrid,ygrid,ycenter,xcenter,sigma)

#read(The unread procar file, the width of gaussians (0=points), ngPnts=number of grid points)
def read(procarFilename,sigma=0.1,ngPnts=1000):
    kpoints=list()
    weights=list()

    procar=open(procarFilename)

    procar.readline()

    [Nkpoints,Nbands,Nions]=map(lambda x:int(x.split(":")[1]),procar.readline().split("#")[1:])
    skbk=procar.tell()
    for i in range(5):
        procar.next()
    Norbs=len(procar.next().split())-1
    procar.seek(skbk)

    #Parsy Parsy noitch.
    print "Parsing PROCAR, this may take a moment..."
    occupancies=zeros([Nkpoints,Nbands,Norbs])
    energies = zeros([Nkpoints,Nbands])
    for i in range(Nkpoints):
        print "On Kpoint %d of %d"%(i+1,Nkpoints)
        k,w=grabKP(procar)
        for j in range(Nbands):
            e=grabE(procar)
            o,orbLabels=grabOatB(procar,Nions)
            occupancies[i,j] = o
            energies[i,j] = e
        kpoints.append(k)
        weights.append(w)

    #Gaussian summations
    eGrid=linspace(energies.min()-1.0,energies.max()+1,ngPnts)
    ocGrids=zeros([Norbs,ngPnts]) #occupancy grid to sum the band-gaussians over

    if sigma>0:
        #Sum up gaussians
        for kp in range(Nkpoints): 
            for b in range(Nbands):
                for o in range(Norbs):
                    ocGrids[o] = sumgauss( energies[kp,b] , occupancies[kp,b,o]*weights[kp]/Nkpoints ,eGrid,ocGrids[o],sigma)
    else:
        #Sump up points
        for kp in range(Nkpoints): 
            for b in range(Nbands):
                for o in range(Norbs):
                    a=digitize([energies[kp,b]],eGrid)
                    ocGrids[o,a] += occupancies[kp,b,o]*weights[kp]/Nkpoints

    #lengths
    #      Norbs    ,Nkpoints,Nkpoints,Nkpoints*NBands,Nkpoints*Norbs*Nbands,ngPnts,Norbs*ngPnts
    return orbLabels,kpoints ,weights ,energies       ,occupancies          ,eGrid ,ocGrids

#Parses the PROCAR averaging over Kpoints and Bands
def read_per_atom(procarFilename)
