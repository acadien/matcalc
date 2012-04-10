from scipy import *
#mine
from geometry import cutPeriodicBounds

def readQVo(qvodata,bounds=None):
    q0=qvodata[0].strip()
    if q0!="3":
        print "This script only works with 3 dimensions, instead %s found in QVo file."%q0
    nverts,natoms,one=map(int,qvodata[1].split())
    qvodata=qvodata[2:]
    #if bounds==None:
    verts=array([map(float,i.split()) for i in qvodata[:nverts]])
    #else:
        #verts=array([map(float,i.split()) for i in qvodata[:nverts]])
    #    verts=array([cutPeriodicBounds(array(map(float,i.split())),bounds) for i in qvodata[:nverts]])
        #print verts
        #verts=[i for i in verts if i!=None] #remove the none's           
    polyverts=[map(int,i.split()[1:]) for i in qvodata[nverts:nverts+natoms]]
    polyverts=[[verts[j] for j in i] for i in polyverts]
    #polyverts=[i for i in polyverts if len(i)>0]
    #print polyverts
    return polyverts

def readQVFi(qvfidata,full=True):
    data=[i.split()[1:] for i in qvfidata[1:]]
    rawplanes=[[int(line[0]),int(line[1]),map(float,line[2:])] for line in data]
    nplanes=max(max(zip(*rawplanes)[0]),max(zip(*rawplanes)[1]))+1
    polyplanes=[list() for i in range(nplanes)]
    neighbors=[list() for i in range(nplanes)]
    for [i,j,coefs] in rawplanes:
        polyplanes[i].append(coefs)
        neighbors[i].append(j)
        if full:
            polyplanes[j].append(coefs)
            neighbors[j].append(i)
    return array(polyplanes),neighbors
