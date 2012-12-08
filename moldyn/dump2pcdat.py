#!/usr/bin/python
import sys
from math import *
from scipy import array
#mine
from lammpsIO import dumpReadNext
from rdf import rdf_periodic
from plotstruct import plot_atoms

#Calculates the RDF for a set of atoms in an XDATCAR file and writes it to a file

def usage():
    print "Usage: %s <lammps-dump> <cutoff=10.0> <nbins=1000> <pcfile=\"PCDAT_fromLAMMPS\">"%sys.argv[0].split("/")[-1]

cutoff=10.0
nbins=1000
pcfile="PCDAT_fromLAMMPs"
if len(sys.argv)<2:
    usage()
    exit(0)
if len(sys.argv)>=3:
    cutoff=float(sys.argv[2])
if len(sys.argv)>=4:
    nbins=int(sys.argv[3])    
if len(sys.argv)==5:
    pcfile=sys.argv[4]

dumpfile=sys.argv[1]
dump=open(dumpfile,"r").readlines()

delr=float(cutoff)/nbins
pcors=list()

i=0
cnt=0
while True:
    try:
        dump,bounds,types,atoms,head = dumpReadNext(dump)
        cnt+=1
    except Exception as inst:
        #last of the configurations
        break
    pcors.append(rdf_periodic(array(atoms),array(bounds),cutoff=cutoff,nbins=nbins)[1])
    print cnt
print "Found %d configurations."%cnt
print "Writing RDF of these configurations to %s"%pcfile
head="\nCreated by %s, from %s\n"%(sys.argv[0].split("/")[-1],dumpfile)+"\n"*4+str(nbins)+"\n"*2+str(delr/1E10)+"\n"*4
pcdat=open(pcfile,"w")
pcdat.write(head)
for pcor in pcors:
    dat="\n".join(map(str,pcor))
    pcdat.write(dat+"\n\n")
    
