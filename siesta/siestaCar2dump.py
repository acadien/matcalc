#!/usr/bin/python

import sys
#mine
import utils
import lammpsIO

utils.usage(["<siesta MD_CAR>","<optional: dumpfile, default:prefix.dat>"],1,2)

inputFile = sys.argv[1]

if len(sys.argv)==3:
    outputFile = sys.argv[2]
else:
    outputFile = ".".join(inputFile.split(".")[:-1]+["dat"])

mdcar = open(inputFile,"r")
opdat = open(outputFile,"w")

count = 0
three = range(3)
while True:
    mdcar.readline() #header
    mdcar.readline() #volume ratio
    v1 = map(float,mdcar.readline().split())
    v2 = map(float,mdcar.readline().split())
    v3 = map(float,mdcar.readline().split())

    nAtom = int(mdcar.readline())
    mdcar.readline() #direct
    atoms = [map(float,mdcar.readline().split()) for i in range(nAtom)]
    if len(atoms[nAtom-1])!=3:
        break
    types = [0 for i in range(nAtom)]

    xhi,yhi,zhi,xy,xz,yz = lammpsIO.basis2lohi([v1,v2,v3])
    bf=[xhi,yhi,zhi]
    atoms = [map(lambda x: (a[x]+0.5)*bf[x],three) for a in atoms]
    
    head= "ITEM: TIMESTEP\n%d\n"%count
    head+="ITEM: NUMBER OF ATOMS\n%d\n"%nAtom
    head+="ITEM: BOX BOUNDS\n"
    head+=" 0.0000  % 6.6f\n"%(xhi)
    head+=" 0.0000  % 6.6f\n"%(yhi)
    head+=" 0.0000  % 6.6f\n"%(zhi)
    head+="ITEM: ATOMS id type xs ys zs\n"
    opdat.write(head)

    atoms = zip(*atoms)
    opdat.writelines([" ".join(map(str,line))+"\n" for line in zip(range(nAtom),types,atoms[0],atoms[1],atoms[2])])
        
    count += 1
print count
