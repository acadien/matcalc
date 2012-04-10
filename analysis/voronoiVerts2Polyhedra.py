#!/usr/bin/python

import sys
#mine
from voronoiIO import *
from poscarIO import readposcar
from struct_tools import points2polyhedron

def usage():
    print "Usage: %s <qvoronoi o-file (Qo_output)> <qvoronoi Fi-file (QFi_output> <POSCAR>"%sys.argv[0].split("/")[-1]

"""
Output is a file with this format:
<N> //Number of Polyhedron
<M0>//Number of Polygons for Polyhedron i
v0_0 v0_1 v0_2 ... v0_L0
v1_0 v1_1 ... v1_L1
...
vM0_0 vM0_1 vM0_2 vM0_3 ... vM0_LM0
<M1>
v0_0 ...
...
vM1_0 ...
...
...
<MN>
v0_0 ...
...
vMN_0 ...
<Nv> //Number of Vertices
<v0xyz>
<v1xyz>
<v2xyz>
"""

if len(sys.argv)<4:
    usage()
    exit(0)

qvodata=open(sys.argv[1],"r").readlines()
qvfidata=open(sys.argv[2],"r").readlines()
poscar=open(sys.argv[3],"r").readlines()

[v1,v2,v3,atypes,ax,ay,az,head,poscar] = readposcar(poscar)
basis=array([v1,v2,v3])
atoms=zip(ax,ay,az)
bounds=array([basis[0][0],basis[1][1],basis[2][2]])

polyverts=readQVo(qvodata)#,bounds)
polyplanes,neighbors=readQVFi(qvfidata)

#Chop off polyhedra of points outside simulation
polyplanes=polyplanes[:len(atoms)]

polyhedra=[points2polyhedron(verts,planes,plotting=False) for verts,planes in zip(polyverts,polyplanes)]

polyverts=readQVo(qvodata)
polyplanes,neighbors=readQVFi(qvfidata)

polyhedra=[points2polyhedron(verts,planes,plotting=False) for verts,planes in zip(polyverts[0:100],polyplanes[0:100])]

data=[str(len(polyhedra))]#Number of polyhedra, should be # of atoms
verts=list(set([i for i in polyhedra]))
for polyhedron in polyhedra:
    data.append(str(len(polyhedron)))#Number of polygons
    for polygon in polyhedron:
        data.append(" ".join(map(str,polygon)))

open("polyhedra.dat","w").writelines([i+"\n" for i in data])
    
