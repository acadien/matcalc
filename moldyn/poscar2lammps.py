#!/usr/bin/python

import sys
from poscarIO import readposcar
from numpy import *
#mine
from struct_tools import mag,ang,flatten

def usage():
    print "%s <in:POSCAR-file> <out:LAMMPS config-file>"%(sys.argv[0].split("/")[-1])

if len(sys.argv)<3:
    usage()
    exit(0)

try:
    poscar=open(sys.argv[1],"r").readlines()
    lammps=open(sys.argv[2],"w")
except IOError:
    print "Error: Unable to open files."
    usage()
    exit(0)

[v1,v2,v3,atypes,ax,ay,az,head,poscar] = readposcar(poscar,frac_coord=True)
v1,v2,v3=map(array,[v1,v2,v3])
atoms=zip(ax,ay,az)
#Convert from POSCAR style basis vectors to LAMMPS style boundaries.
mv1,mv2,mv3=map(mag,[v1,v2,v3]) #magnitude

origin=array([0,0,0])
A=ang(v2,origin,v3)
B=ang(v1,origin,v3)
C=ang(v1,origin,v2)

#Only 14 digits of accuracy.
xhi= mv1
xy= mv2*cos(C)
xz= mv3*cos(B)
yhi=(mv2**2-xy**2)**0.5
yz= (mv2*mv3*cos(A)-xy*xz)/yhi
zhi=(mv3**2-xz**2-yz**2)**0.5

basis=matrix([[xhi,0,0],[xy,yhi,0],[xz,yz,zhi]])

N=sum(atypes)
Ntypes=len(atypes)
types=[i for i in flatten([[i]*v for i,v in enumerate(atypes)])]
data="#Generated by %s from %s\n\n"%(sys.argv[0].split("/")[-1],sys.argv[1])
data+="%d atoms\n"%N
data+="%d atom types\n\n"%Ntypes
data+="0.0000  %4f  xlo xhi\n"%xhi
data+="0.0000  %4f  ylo yhi\n"%yhi
data+="0.0000  %4f  zlo zhi\n"%zhi
data+="% 4f % 4f % 4f xy xz yz\n\n"%(xy,xz,yz)
data+="Atoms\n\n"
for i,(t,atom) in enumerate(zip(types,atoms)):
    x,y,z=(atom*basis).tolist()[0]
    data+="\t %d \t %d \t % 5f\t% 5f\t% 5f\n"%(i+1,t+1,x,y,z)
lammps.write(data)
lammps.close()
