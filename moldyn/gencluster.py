#!/usr/bin/python

import sys
import random
from numpy import matrix

minr=1.6

if len(sys.argv) < 8:
    print "Usage: %s <Natoms> <genx low,high> <geny low,high> <genz low,high> <box boundx> <box boundy> <box boundz>" % (sys.argv[0])
    print "Generates N atoms in random locations within the 'generate' bounding box given, and not with in a radius of %g from any other atoms generated. Than inverts the atomic locations over the bounding box to obtain the fractional coordinates." % (minr)
    print "Example: %s 3,4 1.5,3.5 1.5,3.5 1.5,3.5 5 5 5" % (sys.argv[0])
    exit(0)

Natoms=map(int,sys.argv[1].split(","))
genx=map(float,sys.argv[2].split(","))
geny=map(float,sys.argv[3].split(","))
genz=map(float,sys.argv[4].split(","))
bndx=float(sys.argv[5])
bndy=float(sys.argv[6])
bndz=float(sys.argv[7])

xs=list()
ys=list()
zs=list()

for i in range(sum(Natoms)):
    cnt=0
    fnd=False
    while fnd!=True and cnt<1E6:
        cnt+=1
        x=random.random()*(genx[1]-genx[0])+genx[0]
        y=random.random()*(geny[1]-geny[0])+geny[0]
        z=random.random()*(genz[1]-genz[0])+genz[0]
        fnd=True
        for a,b,c in zip(xs,ys,zs):
            if ((x-a)**2+(y-b)**2+(z-c)**2)**0.5 < minr:
                fnd=False
                break
    if fnd==False:
        print "Unable to make a random arrangement of atoms within the space given, expand the bounds and try again."
        exit(0)
    xs.append(x)
    ys.append(y)
    zs.append(z)
    print "% 7.7g % 7.7g % 7.7g" % (x/bndx,y/bndy,z/bndz)

