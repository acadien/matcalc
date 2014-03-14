#!/usr/bin/python
import scipy
import subprocess
import pylab as pl
import numpy as np
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
    occ = np.array(map(float,procar.readline().split()[1:]))
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
    occupancies=np.zeros([Nkpoints,Nbands,Norbs])
    energies = np.zeros([Nkpoints,Nbands])
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
    eGrid=np.linspace(energies.min()-1.0,energies.max()+1,ngPnts)
    ocGrids=np.zeros([Norbs,ngPnts]) #occupancy grid to sum the band-gaussians over

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
                    a=np.digitize([energies[kp,b]],eGrid)
                    ocGrids[o,a] += occupancies[kp,b,o]*weights[kp]/Nkpoints

    #lengths
    #      Norbs    ,Nkpoints,Nkpoints,Nkpoints*NBands,Nkpoints*Norbs*Nbands,ngPnts,Norbs*ngPnts
    return orbLabels,kpoints ,weights ,energies       ,occupancies          ,eGrid ,ocGrids

#Parses the PROCAR averaging over Kpoints and Bands, giving the local (or Site-Projected) Density of States. LDOS
def readLDOS(procarFilename,sigma=0.1,nGridPoints=1000):
    
    #Parse the header
    pcarF = open(procarFilename,"r")
    pcarF.next()
    nKpoint,nBand,nIon=map(lambda x:int(x.split(":")[1]),pcarF.next().split("#")[1:])
    pcarF.close()

    #Find the starting spots for the band listings
    locBandOccData= subprocess.check_output("grep -b band.*energ.*occ.* %s"%procarFilename,shell=True).split("\n")[:-1]
    locs,bandOccData=zip(*[i.split(":") for i in locBandOccData])
    locs=map(int,locs)
    bandEnergies,occTot=zip(*[map(float,[i.split()[4],i.split()[-1]]) for i in bandOccData])
    bandEnergies=np.asarray(bandEnergies)

    #Find the occupancies for each Ion
    pcarF = open(procarFilename,"r")
    occupancies=np.zeros([nIon,nBand])
    for i in range(len(locs)):
        pcarF.seek(locs[i])
        pcarF.next()
        pcarF.next()
        pcarF.next()
        occupancies[:,i%nBand] += np.asarray([float(pcarF.next().split()[-1]) for j in range(nIon)])
    occupancies/=nKpoint

    #Setup the band energies grid
    mx=max(bandEnergies)+1.0
    mn=min(bandEnergies)-1.0
    bandGrid=np.asarray([float(i)/nGridPoints*(mx-mn)+mn for i in range(nGridPoints)])

    #Setup the occupancy grid and fill with gaussians
    occGrid=np.zeros([nIon,nGridPoints])
    for i in range(nIon):
        #Sums up gaussians accross occGrid on each ion
        [sumgauss(bandEnergies[j],occupancies[i,j],bandGrid,occGrid[i,:],sigma) for j in range(nBand)]
    
    #      float,float      nGP      nIon*nGP
    return nIon,nGridPoints,bandGrid,occGrid
