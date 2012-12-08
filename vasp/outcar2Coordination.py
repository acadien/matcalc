#!/usr/bin/python

import sys
from math import *
import pylab as pl
from scipy import *
import numpy as np
#mine
from geometry import atomsAtLength
from voronoiNeighbors import voronoiNeighbors
from struct_tools import *
import coordinationIO

def outcarCoordination(outcarfile,maxls,voronoiEnable):
    outcar=open(outcarfile,"r")

    #Grab ion types
    while True:
        line=outcar.readline()
        if "ions per type" in line:
            break
    atypes=map(int,line.split("=")[1].split())

    tCoordNumbers=zeros([len(maxls),sum(atypes)])

    #Grab basis vectors
    while True:
        line=outcar.readline()
        if "direct lattice vectors" in line:
            break
    basis=array([map(float,outcar.readline().split()[:3]) for i in range(3)])
    lengths=array([basis[0][0],basis[1][1],basis[2][2]])

    absMax=50
    if voronoiEnable:
        avgCNhist=zeros([len(maxls)+1,absMax])
    else:
        avgCNhist=zeros([len(maxls),absMax])

    #Grab atom positions and perform adf calculation
    count=0
    posit=False
    for line in outcar:
        if posit:
            if "--" in line:
                if len(atoms)==0:
                    continue
                else:
                    #Analysis
                    if max(zip(*atoms)[0])<1.01:
                        atoms=dot(array(atoms),basis)
                    else:
                        atoms=array(atoms)

                    atLFullNeighbors=[neighbors(atoms,array([[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]),maxlen,style='full') for maxlen in maxls]
                    atLCoordNumbers=map(lambda x:map(len,x),atLFullNeighbors)

                    if voronoiEnable:
                        vFullNeighbors=voronoiNeighbors(atoms=atoms,basis=basis,atypes=atypes,style='full')
                        vFullNeighbors=[[j for j in vFullNeighbors[i] if dist_periodic(atoms[i],atoms[j],lengths)<voronCap] for i in range(len(vFullNeighbors))]
                        vCoordNumbers=map(len,vFullNeighbors)
                        atLCoordNumbers.append(vCoordNumbers)
                    
                    hist=[[i.count(j) for j in range(absMax)] for i in atLCoordNumbers]
                    avgCNhist+=array(hist)
                    
                    print count
                    posit=False
            else:
                atoms.append(map(float,line.split()[:3]))
        elif "POSITION" in line:
            atoms=list()
            posit=True
            count+=1
    
    avgCNhist/=count
    return avgCNhist.tolist()


if __name__=="__main__":

    if len(sys.argv)<3:
        print "Usage:"
        print sys.argv[0]+" <outcar> <maxlen0> <maxlen1> <Voronoi,Cap>..."
        exit(0)

    outcarfile=sys.argv[1]
    
    maxls=list()
    voronoiEnable=False
    voronCap=0
    for i in sys.argv[2:]:
        try:
            maxls.append(float(i))
        except ValueError:
            if i[0] in ["v","V"]:
                voronoiEnable=True
                voronCap=float(i.split(",")[1])
            else:
                raise ValueError

    atLCNhist=outcarCoordination(outcarfile,maxls,voronoiEnable)

    colors=["red","blue","green","purple","yellow"]

    if voronoiEnable:
        labels=["<"+str(l) for l in maxls]+["Voronoi <%g"%voronCap]
    else:
        labels=[str(l) for l in maxls]

    header=""
    avgs=list()
    mxcn=0
    mncn=50
    for k,hist in enumerate(atLCNhist):
        avgs.append(sum([i*v for i,v in enumerate(hist)])/sum(hist))
        header+="For bonds of type %s : CN=%g.  "%(labels[k],avgs[-1])
        hist=array(hist)
        a=[i for i in np.where(hist != 0)][0]
        mxcn=max([mxcn,max(a)])
        mncn=min([mncn,min(a)])
    mxcn+=4
    mncn=max([0,mncn-3])
    #Write this information to a file
    try:
        num=sys.argv[1].split("_")[1]
    except IndexError:
        num=None

    if num:
        cnfile="coordination_%s%s.data"%(num,"_"+"_".join(sys.argv[2:]))
    else:
        cnfile="coordination%s.data"%("_"+"_".join(sys.argv[2:]))

    print header
    print "Writing %s."%cnfile

    coordinationIO.write(cnfile,header,0,mxcn,labels,avgs,atLCNhist)
