#!/usr/bin/python
import sys

#This script should be used to copy as POSCAR structure multiple times, in effect growing the system
#There should be no newlines in the input POSCAR file
#The new POSCAR is written to stdout

if not(len(sys.argv) in [5,6]):
    print "Usage:"
    print "%s <dir/POSCAR file> <copyx> <copyy> <copyz> <1:ploton>"%sys.argv[0]
    exit(0)

fil=open(sys.argv[1],"r")
copyx=int(sys.argv[2])
copyy=int(sys.argv[3])
copyz=int(sys.argv[4])

ploton=0
if len(sys.argv)==6:
    ploton=1
    from enthought.mayavi.mlab import figure,points3d,show

xs=list()
ys=list()
zs=list()

#Capture the header, before the atomic positions
head="This POSCAR manipulated by %s script.  "%(sys.argv[0])
head+=fil.readline()
head+=fil.readline()

#Begin by rescaling the bounding vectors by the copy parameters
xvec=[float(i) for i in fil.readline().split()]
yvec=[float(i) for i in fil.readline().split()]
zvec=[float(i) for i in fil.readline().split()]
head+=" ".join([str(i*copyx) for i in xvec])+"\n"
head+=" ".join([str(i*copyy) for i in yvec])+"\n"
head+=" ".join([str(i*copyz) for i in zvec])+"\n"

#Get the number of each atom type
atomtypes=[int(i) for i in fil.readline().split()]
head+=" ".join([str(i*copyx*copyy*copyz) for i in atomtypes])+"\n"

#Check for direct scaling
if fil.readline() in ["Direct","direct"]:
 print "Error: This script only works with Direct scaling, not cartesian."
 exit(0)
head+="Direct"

#Grab the atom locations
for at in atomtypes:
    xs.append(list())
    ys.append(list())
    zs.append(list())
    for n in range(at):
        [x,y,z]=[float(i) for i in fil.readline().split()]
        xs[-1].append(x)
        ys[-1].append(y)
        zs[-1].append(z)

if ploton==1:
    figure()
    for i in range(len(atomtypes)):
        points3d([xs[i][j]*xvec[0]+ys[i][j]*yvec[0]+zs[i][j]*zvec[0] for j in range(len(xs[i]))],[xs[i][j]*xvec[1]+ys[i][j]*yvec[1]+zs[i][j]*zvec[1] for j in range(len(ys[i]))],[xs[i][j]*xvec[2]+ys[i][j]*yvec[2]+zs[i][j]*zvec[2] for j in range(len(zs[i]))])

#Shift the atoms
for ix in range(copyx):
    for iy in range(copyy):
        for iz in range(copyz):
            if ix==iy==iz==0:
                continue
            for i in range(len(atomtypes)):
                for j in range(atomtypes[i]):
                    xs[i].append(xs[i][j]+ix)
                    ys[i].append(ys[i][j]+iy)
                    zs[i].append(zs[i][j]+iz)
                
for i in range(len(atomtypes)):
    xs[i]=[j/copyx for j in xs[i]]
    ys[i]=[j/copyy for j in ys[i]]
    zs[i]=[j/copyz for j in zs[i]]

print head

for i in range(len(atomtypes)):
    for a,b,c in zip(xs[i],ys[i],zs[i]):
        print "  %5.5f %5.5f %5.5f"%(a,b,c)

for i in range(3):
    xvec[i]*=copyx
    yvec[i]*=copyy
    zvec[i]*=copyz

if ploton==1:
    figure()
    for i in range(len(atomtypes)):
        points3d([xs[i][j]*xvec[0]+ys[i][j]*yvec[0]+zs[i][j]*zvec[0] for j in range(len(xs[i]))],[xs[i][j]*xvec[1]+ys[i][j]*yvec[1]+zs[i][j]*zvec[1] for j in range(len(ys[i]))],[xs[i][j]*xvec[2]+ys[i][j]*yvec[2]+zs[i][j]*zvec[2] for j in range(len(zs[i]))])
    show()
