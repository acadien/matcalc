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

def cutPeriodicBounds(pnt,bnds):
    code = """
    int tot=0;
    for(int i=0;i<3;i++){
        if(pnt[i]>=bnds[i]){
            pnt[i]=bnds[i];
            tot++;
        }
        if(pnt[i]<0){
            pnt[i]=0;
            tot++;
        }
    }
    if(tot==3)
        return_val=-1;
    else
        return_val=0;
    """
    if weave.inline(code,['pnt','bnds']) == -1:
        return None
    return pnt

def applyPeriodicBounds(pnt,bnds):
    code = """
    for(int i=0;i<3;i++){
        while(pnt[i]>=bnds[i])
            pnt[i]-=bnds[i];
        while(pnt[i]<0)
            pnt[i]+=bnds[i];
    }
    """
    weave.inline(code,['pnt','bnds'])      
    return pnt
