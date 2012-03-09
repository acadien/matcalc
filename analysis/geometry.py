#!/usr/bin/python

from scipy import *
#mine
from struct_tools import *

#Given a list of atoms and periodic boundary conditions, returns a list of pairs of indeces that are atoms at that distance
#key part: fabs(dist(a,b)-length) < err => true
def atomsAtLength(atoms,halfNeighbors,length,err=0.05,periodic=False,bounds=None):
    keypairs=list()
    for iNeighb,Neighbs in enumerate(halfNeighbors):
        atomi=atoms[iNeighb]
        for jNeighb in Neighbs:
            atomj=atoms[jNeighb]
            if periodic:
                if fabs(dist_periodic(atomi,atomj,bounds)-length) < err:
                    keypairs.append([iNeighb,jNeighb])
            else:
                if fabs(dist(atomi,atomj)-length) < err:
                    keypairs.append([iNeighb,jNeighb])

    return keypairs
