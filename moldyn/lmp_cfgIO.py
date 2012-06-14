#!/usr/bin/python

import re
from scipy import array

#Reads in an atomic configuration and returns relavent bits.
def readcfg(cfgfile):

    nAtomsItem=re.compile(".*NUMBER OF ATOMS.*")
    atomsItem=re.compile(".*ATOMS id type.*")
    boxBoundItem=re.compile(".*BOX BOUNDS.*")

    for i,line in enumerate(cfgfile):
        if "ITEM:" in line:
            if nAtomsItem.match(line):
                inatom=i+1
            elif atomsItem.match(line):
                iatom=i+1
            elif boxBoundItem.match(line):
                ibb=i+1

    Natom=int(cfgfile[inatom])

    #bounds
    bounds=[map(float,i.split()) for i in cfgfile[ibb:ibb+3]]
    
    #type,positions,velocities
    [types,ax,ay,az,avx,avy,avz]=zip(*[map(float,line.split()[1:]) for line in cfgfile[iatom:iatom+Natom]])
    types=map(int,types)

    axyz=array(zip(ax,ay,az))
    avxyz=array(zip(avx,avy,avz))

    return bounds,axyz,avxyz,types
