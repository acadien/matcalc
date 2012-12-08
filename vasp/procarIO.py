#!/usr/bin/python
import scipy
from scipy import array,zeros

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
    occ=zeros(Norbs)
    occ=map(lambda x:sum(x),zip(*[map(float,procar.readline().split()[1:]) for i in range(Nions)]))
    return occ,header

def read(procarFilename):
    energy=list()
    kpoints=list()
    weights=list()
    occupancies=list()#per kpoint, averaged over sites
    avgoc=list()#occupancy averaged over kpoints

    procar=open(procarFilename)
    procar.readline()
    [Nkpoints,Nbands,Nions]=map(lambda x:int(x.split(":")[1]),procar.readline().split("#")[1:])
    print "Parsing PROCAR, this may take a moment..."
    for i in range(Nkpoints):
        print "On Kpoint %d of %d"%(i+1,Nkpoints)
        k,w=grabKP(procar)
        occupancies.append(list())
        energy.append(list())
        for j in range(Nbands):
            e=grabE(procar)
            o,orbLabels=grabOatB(procar,Nions)
            occupancies[i].append(o)
            energy[i].append(e)
        kpoints.append(k)
        weights.append(w)

    energy=array(energy)
    occupancies=array([array(i).T for i in occupancies]) #occupancies[Nkpoints][Nbands][orbitals]

    avgoc=[array(v)/Nkpoints for i,v in enumerate(occupancies)]#multiply by weights
    avgoc=scipy.sum(avgoc,axis=0)
    #avgoc=avgoc.T
    avgen=map(lambda x:sum(x)/Nkpoints,zip(*energy))

    #      Norbs    ,Nkpoints,Nkpoints,Nkpoints*NBands,Nkpoints*Norbs*Nbands,NBands,Norbs*Nbands
    return orbLabels,kpoints ,weights ,energy         ,occupancies          ,avgen ,avgoc       
