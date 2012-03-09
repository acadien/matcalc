#!/usr/bin/python

import sys
import random

def usage():
    print "Usage: %s <[SC,BCC,FCC]> Ncx Ncy Ncz Cx Cy Cz <filename> [ploton=1 else 0]" % (sys.argv[0])
    print "Generates a bunch of atoms on the lattice requested (SC/BCC/FCC), the cell is duplicated Ncx*Ncy*Ncz times with lattice constants Cx/Cy/Cz."
    print "Example: %s FCC 2 2 2 3.165 3.165 3.165" % (sys.argv[0])
    print "The number of atoms, lattice vectors and atomic positions are printed out."
    print "Dumps the atomic configuration to a file <filename> which is readable by LAMMPS through the read_data command."
    exit(0)

bcc=[[0,0,0],[0.5,0.5,0.5]]#[1.0,0.0,0.0],[0.0,1.0,0.0]]
fcc=[[0,0,0],[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]]
sc=[[0,0,0]]#,[1,0,0],[0,1,0],[0,0,1]]
lats={'SC':sc,'BCC':bcc,'FCC':fcc}

if len(sys.argv) != 10:
    usage()

ltype=sys.argv[1]
if not(ltype in lats.keys()):
    print "Bad input."
    usage()

[Ncx,Ncy,Ncz]=map(int,sys.argv[2:5])
[Cx,Cy,Cz]=map(float,sys.argv[5:8])
ploton=int(sys.argv[9])

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

#Set the bounds
bndx=Ncx*Cx
bndy=Ncy*Cy
bndz=Ncz*Cz

#Center the atoms about 1/2 the bounds
xshft=bndx/2.0-sum(xs)/len(xs)
yshft=bndy/2.0-sum(ys)/len(ys)
zshft=bndz/2.0-sum(zs)/len(zs)
xs=map(lambda x:x+xshft,xs)
ys=map(lambda y:y+yshft,ys)
zs=map(lambda z:z+zshft,zs)

N=len(xs)
print "Lattice has %d atoms.\n" % N

#Plotting
if ploton==1:
    from enthought.mayavi import mlab
    fig=mlab.gcf()
    mlab.points3d(xs,ys,zs)
    mlab.plot3d([0,bndx],[0,0],[0,0])
    mlab.plot3d([0,0],[0,bndy],[0,0])
    mlab.plot3d([0,0],[0,0],[0,bndz])
    mlab.show()

#Ask before overwriting a file.
try:
    val=raw_input("Continue? (y/n): ")
    if not(val in ('y','Y')):
        print "Okey dokey, stopping."
        exit(0)
    ofil=open(sys.argv[8],"w")
except ValueError:
    exit(0)

#Dump the atomic coordinates etc to the file, make it readable for lammps
line  = "#This is an atomic configuration file for a %d atoms in a %s configuration, the atoms are centered on a lattice made for periodic boundary conditions.\n\n" % (N,ltype)
line += "\t%d atoms\n" % (N)
line += "\n"
line += "\t1 atom types\n"
line += "\n"
line += "0.0 %4.4g  xlo xhi\n" % (bndx)
line += "0.0 %4.4g  ylo yhi\n" % (bndy)
line += "0.0 %4.4g  zlo zhi\n" % (bndz)
line += "0.0 0.0 0.0 xy xz yz\n"
line += "\nAtoms\n\n"
for i in range(N):
    line += "\t%d\t%d\t%6.6g\t%6.6g\t%6.6g\n" % (i+1,1,xs[i],ys[i],zs[i])
ofil.write(line)


