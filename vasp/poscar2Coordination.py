#!/usr/bin/python

import sys
from math import *
import pylab as pl
from scipy import *
#mine
import poscarIO
from geometry import atomsAtLength
from voronoiNeighbors import voronoiNeighbors
from struct_tools import *
import coordinationIO

if len(sys.argv)<3:
    print "Usage:"
    print sys.argv[0]+" <poscar> <maxlen0> <maxlen1> <Voronoi,Cap>..."
    exit(0)

poscar = open(sys.argv[1],"r").readlines()
maxls=list()
voronEnable=False
voronCap=0
for i in sys.argv[2:]:
    try:
        maxls.append(float(i))
    except ValueError:
        if i[0] in ["v","V"]:
            voronEnable=True
            voronCap=float(i.split(",")[1])
        else:
            raise ValueError

[basis,atypes,atoms,head,poscar] = poscarIO.read(poscar)
atoms=array(atoms)
lengths=array([basis[0][0],basis[1][1],basis[2][2]])

#vHalfNeighbors=voronoiNeighbors(atoms=atoms,basis=basis,atypes=atypes,style='half')
#atLHalfNeighbors=[neighbors(atoms,array([[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]),maxlen,style='half') for maxlen in maxls]
atLFullNeighbors=[neighbors(atoms,array([[0,basis[0][0]],[0,basis[1][1]],[0,basis[2][2]]]),maxlen,style='full') for maxlen in maxls]
vFullNeighbors=voronoiNeighbors(atoms=atoms,basis=basis,atypes=atypes,style='full')
#atLFullNeighbors=map(half2full,atLHalfNeighbors)
#vFullNeighbors=half2full(vHalfNeighbors)
vFullNeighbors=[[j for j in vFullNeighbors[i] if dist_periodic(atoms[i],atoms[j],lengths)<voronCap] for i in range(len(vFullNeighbors))]

atLCoordNumbers=map(lambda x:map(len,x),atLFullNeighbors)
vCoordNumbers=map(len,vFullNeighbors)

#To include Voronoi (capped at max length) uncomment:
if voronEnable:
    atLCoordNumbers.append(vCoordNumbers)

absMin=0
absMax=50
mncn=max(min(map(min,atLCoordNumbers))-3,absMin)
mxcn=min(max(map(max,atLCoordNumbers))+4,absMax)

colors=["red","blue","green","purple","yellow"]
atLCNhist=list()

if voronEnable:
    labels=["<"+str(l) for l in maxls]+["Voronoi <%g"%voronCap]
else:
    labels=[str(l) for l in maxls]

header=""
avgs=list()
for k,coordNumbers in enumerate(atLCoordNumbers):
    #For the bond length in question, how many bonds of length in question?
    CNhist=zeros(absMax)
    
    for cn in coordNumbers:
        CNhist[cn]+=1
    avgs.append(sum([i*v for i,v in enumerate(CNhist)])/sum(CNhist))
    header+="For bonds of type %s : CN=%g.  "%(labels[k],avgs[-1])
    atLCNhist.append(list(CNhist))
    c=colors[k]
    #pl.bar(range(mncn,mxcn),CNhist[mncn:mxcn],width=1.0,bottom=0,color=c,alpha=0.75,label=labels[k])
    #pl.title(sys.argv[1])
    #pl.plot([avg,avg],[0,max(CNhist)],c='black',lw=2,ls="-")

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

coordinationIO.write(cnfile,header,mncn,mxcn,labels,avgs,atLCNhist)
