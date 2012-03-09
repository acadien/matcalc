#!/usr/bin/python

import sys
import random

if len(sys.argv) != 9:
    print "Usage: %s <[SC,BCC,FCC]> Ncx Ncy Ncz Cx Cy Cz [plot:0,1]" % (sys.argv[0])
    print "Generates a bunch of atoms on the lattice requested (SC/BCC/FCC), the cell is duplicated Ncx*Ncy*Ncz times with lattice constants Cx/Cy/Cz."
    print "Example: %s FCC 2 2 2 3.165 3.165 3.165" % (sys.argv[0])
    print "The number of atoms, lattice vectors and atomic positions are printed out."
    exit(0)

ltype=sys.argv[1]
[Ncx,Ncy,Ncz]=map(int,sys.argv[2:5])
[Cx,Cy,Cz]=map(float,sys.argv[5:8])
ploton=int(sys.argv[8])

bcc=[[0,0,0],[0.5,0.5,0.5]]
fcc=[[0,0,0],[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]]
sc=[[0,0,0]]#,[1,0,0],[0,1,0],[0,0,1]]
lats={'SC':sc,'BCC':bcc,'FCC':fcc}

if ltype=='SC':
    Ncx+=1
    Ncy+=1
    Ncz+=1

lat=lats[ltype]

#Start making atom coordinates
xs=list()
ys=list()
zs=list()
for xi in range(Ncx):
    for yi in range(Ncy):
        for zi in range(Ncz):
            for base in lat:
                xs.append(Cx*(xi+base[0]))
                ys.append(Cy*(yi+base[1]))
                zs.append(Cz*(zi+base[2]))
#Remove any duplicates
uniques=list(set(zip(xs,ys,zs)))
xs=map(lambda x:x[0],uniques)
ys=map(lambda x:x[1],uniques)
zs=map(lambda x:x[2],uniques)

#Calculate the total lattice vectors, requires a buffer region
if ltype=='SC' and Ncx%2==0:
    buff=1.0
else:
    buff=0

bndx=(Ncx+buff)*Cx
bndy=(Ncy+buff)*Cy
bndz=(Ncz+buff)*Cz

#Center the atoms about 1/2 the bounds
xshft=bndx/2.0-sum(xs)/len(xs)
yshft=bndy/2.0-sum(ys)/len(ys)
zshft=bndz/2.0-sum(zs)/len(zs)

xs=map(lambda x:x+xshft,xs)
ys=map(lambda y:y+yshft,ys)
zs=map(lambda z:z+zshft,zs)

print "N=%d" % (len(xs))
print "Lattice vectors:"
print "% 4.4g % 4.4g % 4.4g" % (bndx,0,0)
print "% 4.4g % 4.4g % 4.4g" % (0,bndy,0)
print "% 4.4g % 4.4g % 4.4g" % (0,0,bndz)
print "Atomic locations:"
for i in zip(xs,ys,zs): print "% 6.6g % 6.6g % 6.6g" % (i)

if ploton==1:
    from enthought.mayavi import mlab
    fig=mlab.gcf()
    mlab.points3d(xs,ys,zs)
    mlab.plot3d([0,bndx],[0,0],[0,0])
    mlab.plot3d([0,0],[0,bndy],[0,0])
    mlab.plot3d([0,0],[0,0],[0,bndz])
    mlab.show()

