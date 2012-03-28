#!/usr/bin/python
import sys,subprocess
from numpy import *
import pylab as pl
from scipy import weave
from scipy.weave import converters
#mine
from poscar2qvoronoi import poscar2qvoronoi
from voronoiIO import readQVFi,readQVo
from struct_tools import *

def applyPeriodicBounds(pnt,bnds):
    code = """
    for(int i=0;i<3;i++){
        if(pnt[i]>=bnds[i])
            pnt[i]-=bnds[i];
        if(pnt[i]<0)
            pnt[i]+=bnds[i];
    }
    """
    weave.inline(code,['pnt','bnds'])      
    #for i in range(3):
    #    if pnt[i]>=bnds[i]:
    #        pnt[i]-=bnds[i]
    #    elif pnt[i]<0:
    #        pnt[i]+=bnds[i]
    return pnt

def roundify(atom):
    return map(lambda x:trunc(x*10.01),atom)

def findatom(aval,alist):
    for i,v in enumerate(alist):
        if dist(v,aval)<0.1:
            return i
    raise ValueError("%s not found in list"%str(aval))

def voronoiNeighbors(**kwargs):
    #Possible Args:
    #Option1: poscar,[style(full/half)]
    #Option2: atoms,basis,atypes,[style(full/half)]

    if "style" in kwargs: style=kwargs["style"]
    else:                 style="Full"

    qvfile,basis,qvatoms,nRealAtoms=poscar2qvoronoi(**kwargs)

    #Get hyperplanes & neighbor list from Fi setting of qvoronoi
    p1 = subprocess.Popen(["cat",qvfile],stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['qvoronoi','Fi'],stdin=p1.stdout,stdout=subprocess.PIPE).stdout
    qvfidata=p2.readlines()
    polyplanes,allNeighbors=readQVFi(qvfidata)

    #Apply boundaries to ghost atoms
    bounds=array([basis[0][0],basis[1][1],basis[2][2]])
    qvatoms=array(qvatoms)

    #Shift ghosts back into simulation-box and re-form neighborlist
<<<<<<< HEAD
    boundedAtoms=array(map(lambda x:applyPeriodicBounds(x,bounds),qvatoms))
    #realAtoms=#[list(roundify(atom)) for atom in qvatoms[:nRealAtoms]]
    realAtoms=qvatoms[:nRealAtoms]
    ghost2real=range(nRealAtoms)+[findatom(array(atom),realAtoms) for atom in boundedAtoms[nRealAtoms:]]

=======
    boundedAtoms=map(lambda x:list(roundify(applyPeriodicBounds(x,bounds))),qvatoms)
    #realAtoms=map(list(roundify(atom)),boundedAtoms)
    realAtoms=boundedAtoms[:nRealAtoms]
    ghost2real=range(nRealAtoms)+[realAtoms.index(atom) for atom in boundedAtoms[nRealAtoms:]]
    """except ValueError:
        for atom in boundedAtoms[nRealAtoms:]:
            try:
                realAtoms.index(roundify(atom))
            except ValueError:
                for i in realAtoms:
                    print i
                print atom
                print roundify(atom)
            exit(0)
        exit(0)
    """
>>>>>>> bf58e679a7f844f24dc34b9bd39c1cb31c60e3a5
    #Build the real neighbor list
    
    neighbors=[list() for i in range(nRealAtoms)]
    for i,neighbs in enumerate(allNeighbors[:nRealAtoms]):
        ir=ghost2real[i]
        for jNeighb in neighbs:
            jr=ghost2real[jNeighb]
            if jr not in neighbors[ir]:
                neighbors[ir].append(jr)

    #Build half-neighbor list if necessary
    hneighbors=[list() for i in range(nRealAtoms)]
    if style in ["h","H","half","Half"]:
        for i,neighbs in enumerate(neighbors):
            hneighbors[i]=[j for j in neighbs if i not in hneighbors[j]]
        return hneighbors
    return neighbors


if __name__=="__main__":
    def usage():
        print "Usage: %s <POSCAR>"%(sys.argv[0])

    if len(sys.argv)<2:
        usage()
        exit(0)

    poscar=open(sys.argv[1],"r").readlines()
    neighbors=voronoiNeighbors(poscar=poscar)

    print neighbors


