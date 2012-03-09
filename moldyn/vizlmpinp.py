#!/usr/bin/python

import sys
import pylab as pl
#Plots atoms taken from an input file for LAMMPS

if len(sys.argv)==1:
    print "Usage: %s <input file>"%sys.argv[0]

ifil=open(sys.argv[1],"r")

Natoms=0
Ntypes=0
for line in ifil:
    if "atoms" in line:
        try:
            Natoms=int(line.split()[0])
        except ValueError:
            pass

    if "atom types" in line:
        try:
            Ntypes=int(line.split()[0])
        except ValueError:
            pass

    if line[0:5]=="Atoms":
        break

if Ntypes==0 or Natoms==0:
    print "Error: Unable to read Number of atoms and/or Number of atom types from LAMMPS input file, please verify the file is formatted correctly."
    exit(0)

atoms=[]
for i in range(Ntypes):
    atoms.append(list())
for line in ifil:
    if len(line)<=5: continue
    line=line.split()
    t=int(line[1])-1 #type counting starts at 1 in LAMMPS
    x=float(line[2])
    y=float(line[3])
    z=float(line[4])
    atoms[t].append([x,y,z])

#3d scatter plot
import mpl_toolkits.mplot3d.axes3d as p3
ax = p3.Axes3D(pl.figure())
for i in range(Ntypes):
    ax.scatter3D(*zip(*atoms[i]))
pl.show()
