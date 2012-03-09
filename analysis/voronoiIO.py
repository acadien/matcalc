from scipy import *

def readQVo(qvodata):
    q0=qvodata[0].strip()
    if q0!="3":
        print "This script only works with 3 dimensions, instead %s found in QVo file."%q0
    nverts,natoms,one=map(int,qvodata[1].split())
    qvodata=qvodata[2:]
    verts=array([map(float,i.split()) for i in qvodata[:nverts]])
    polyverts=[map(int,i.split()[1:]) for i in qvodata[nverts:nverts+natoms]]
    polyverts=[[verts[j] for j in i] for i in polyverts]
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
