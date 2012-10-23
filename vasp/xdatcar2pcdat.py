#!/usr/bin/python
import sys
from math import *
from scipy import array
import pylab as pl
#mine
from paircor import paircor_periodic
from plotstruct import plot_atoms
from duplicate import duplicate26

#Calculates the pair-correlation function for a set of atoms in an XDATCAR file and writes it to a file

def usage():
    print "Usage: %s <xdatcar> <cutoff (A)> <numbins> <opt:pcdat file name>"%sys.argv[0].split("/")[-1]

if len(sys.argv)<4:
    usage()
    exit(0)
if len(sys.argv)==5:
    pcfile=sys.argv[4]
else:
    pcfile="PCDAT_fromX"

xdatfile=sys.argv[1]
xdatcar=open(xdatfile,"r").readlines()
Natoms=int(xdatcar[0].split()[0])
lengths=array(map(lambda x:float(x)*1E10,xdatcar[1].split()[1:4]))
basis=[[1,0,0],[0,1,0],[0,0,1]]
xdatcar=xdatcar[5:]
cutoff=float(sys.argv[2])
numbins=int(sys.argv[3])
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
    pcors.append(paircor_periodic(array(atoms),array(lengths),cutoff=cutoff,nbins=numbins)[1])



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
    
