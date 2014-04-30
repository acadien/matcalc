#!/usr/bin/python

import sys
import outcarIO,poscarIO

#parses an outcar, takes 3 nearest poscar information and uses finitie difference method to calculate velocities. Stick velocities at the end. Can take the place of a poscar to initialize velocities.

def usage():
    print "%s OUTCAR CONFIG"%sys.argv[0].split("/")[-1]

if len(sys.argv)!=3:
    usage()
    exit(0)

outcarFile = sys.argv[1]
configRequest = int(sys.argv[2])

nSteps = outcarIO.nSteps(outcarFile)

#unit is femto, but in 10E-10*10E-15 seconds to account for angstroms in distance
timestep = outcarIO.timestep(outcarFile) * 10E-1

#Error check requested configuration
if configRequest == -1:
    configRequest = nSteps-1
if configRequest > nSteps:
    print "%d configurations available"%nSteps
    exit(0)
if configRequest < 5:
    print "Need more than 5 steps for backwards - finite difference method."
    exit(0)

#Gather poscar information for last set of atoms
basis,atypes,atoms,head,poscar = outcarIO.outcar2poscar(outcarFile,configRequest)

#Convert from fractional units to cartesian
def u2a(u):
    a=list()
    a.append(u[0]*basis[0][0]+u[0]*basis[1][0]+u[0]*basis[2][0])
    a.append(u[0]*basis[0][1]+u[1]*basis[1][1]+u[1]*basis[2][1])
    a.append(u[0]*basis[0][2]+u[2]*basis[1][2]+u[2]*basis[2][2])
    return a

stepb0 = map(u2a,atoms)
stepb1 = map(u2a,outcarIO.outcar2poscar(outcarFile,configRequest-1)[2])
stepb2 = map(u2a,outcarIO.outcar2poscar(outcarFile,configRequest-2)[2])
stepb3 = map(u2a,outcarIO.outcar2poscar(outcarFile,configRequest-3)[2])
stepb4 = map(u2a,outcarIO.outcar2poscar(outcarFile,configRequest-4)[2])

velo5=list()
for b0,b1,b2,b3,b4 in zip(stepb0,stepb1,stepb2,stepb3,stepb4):
    velo5.append( [(25*b0[i]-48*b1[i]+36*b2[i]-16*b3[i]+3*b4[i])/(12.0*timestep) for i in range(3)] )

poscarIO.write("VDATCAR",basis,atoms,atypes,head,velocities=velo5)

