#!/usr/bin/python
import sys
from math import *
from scipy import array
import pylab as pl
#mine
from poscarIO import readposcar
from paircor import paircor_periodic
from plotstruct import plot_atoms

#Calculates the pair-correlation function for a set of atoms in an XDATCAR file and writes it to a file

def usage():
    print "Usage: %s <xdatcar> <POSCAR/CONTCAR (for lattice vectors)> <cutoff (A)> <numbins> <opt:pcdat file name>"%sys.argv[0].split("/")[-1]

cutoff=10.0
numbins=1000
pcfile="PCDAT_fromX"
if len(sys.argv)<3:
    usage()
    exit(0)
if len(sys.argv)>=4:
    cutoff=float(sys.argv[3])
if len(sys.argv)>=5:
    numbins=int(sys.argv[4])    
if len(sys.argv)==6:
    pcfile=sys.argv[5]


xdatfile=sys.argv[1]
pcarfile=sys.argv[2]
v1,v2,v3,d1,d2,d3,d4,d5,d6=readposcar(open(pcarfile,"r").readlines())
basis=[v1,v2,v3]
xdatcar=open(xdatfile,"r").readlines()
Natoms=int(xdatcar[0].split()[0])
lengths=array(map(lambda x:float(x)*1E10,xdatcar[1].split()[1:4]))
xdatcar=xdatcar[5:]

delr=float(cutoff)/numbins
pcors=list()

i=0
cnt=0
atomconfigs=list()
while i < len(xdatcar):
    i+=1
    cnt+=1
    atoms=map(lambda x:map(float,x.split()),xdatcar[i:i+Natoms])
    i+=Natoms
    pcors.append(paircor_periodic(array(atoms),basis,cutoff=cutoff,nbins=numbins)[1])

print "="*50
print "="*50
print "%d Configurations found."%len(pcors)
print "="*50
print "="*50
print "Writing pair distribution of these configurations to %s"%pcfile
head="\nCreated by %s, from %s\n"%(sys.argv[0].split("/")[-1],xdatfile)+"\n"*4+str(numbins)+"\n"*2+str(delr/1E10)+"\n"*4
pcdat=open(pcfile,"w")
pcdat.write(head)
for pcor in pcors:
    dat="\n".join(map(str,pcor))
    pcdat.write(dat+"\n\n")
    
