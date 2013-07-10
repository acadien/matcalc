#!/usr/bin/python
import sys,os

#mine
from outcarIO import outcarReadConfig

#This script writes an OUTCAR configuration to a force database for operating on by potfit

def usage():
    print "Usage:"
    print sys.argv[0].split("/")[-1]+" <force database> <OUTCAR file> <configuration #> <elements (comma separated)> <optional:weight sweight> <optional: scale enable 0>"
    exit(0)

if not(len(sys.argv) in [5,6,7,8]):
    usage()

ocfile = sys.argv[2] #outcar directory

forcefil= sys.argv[1]
fordb = open(forcefil,"r")
grabconfig=int(sys.argv[3])
elems = " ".join(sys.argv[4].split(","))
weight=-1
sweight=-1
scaleEnable=True
if len(sys.argv)>=6:
    weight=int(sys.argv[5])
if len(sys.argv)==7:
    sweight=int(sys.argv[6])
if len(sys.argv)==8:
    if sys.argv[7]==0:
        scaleEnable=False

#get the force configuration number (the total number of configs in the forcefil)
dbcnfgcnt=0 #Force Data Base Configuration Counter
for line in fordb:
    if "#N" in line[0:2]:
        dbcnfgcnt+=1

fordb.close()
fordb = open(forcefil,"a")

TE,stress,basis,atoms,forces,types=outcarReadConfig(ocfile,grabconfig)
#Change atoms to be in cartesian coords instead of fractional
#bT=basis.T
atoms=[basis.dot(atom) for atom in atoms]
natom=len(atoms)

#Number of atoms, use force, header
if scaleEnable:
    line="#N\t%d 1 ifconf=%d Taken From:%s  Config:#%d\n"%(natom,dbcnfgcnt,os.getcwd()+"/"+ocfile,grabconfig)
else:
    line="#N\t%d 1 ifconf=%d Taken From:%s  Config:#%d scaleEnable0\n"%(natom,dbcnfgcnt,os.getcwd()+"/"+ocfile,grabconfig)

#atom types
line += "#C\t"+elems+"\n"

#Bounding box (Angstroms)
line += "#X\t %12.8f  %12.8f %12.8f\n"%tuple(basis[0])
line += "#Y\t %12.8f  %12.8f %12.8f\n"%tuple(basis[1])
line += "#Z\t %12.8f  %12.8f %12.8f\n"%tuple(basis[2])

#Cohesive energy per atom (eV)
line += "#E\t %12.8f\n"%(TE/natom)

#Weight
if weight>0 and sweight<0:
    line += "#W\t %d\n"%weight
elif weight>0 and sweight>0:
    line += "#W\t %d %d\n"%(weight,sweight)

#Stress Tensor
line += "#S\t %12.8f  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f\n"%(stress[0],stress[1],stress[2],stress[3],stress[4],stress[5]) 
            
#Positions
line += "#F\t\n"      
for i in range(natom):
    line += "  %d\t %12.8f  %12.8f  %12.8f   %12.8f  %12.8f  %12.8f\n"%(types[i],atoms[i][0],atoms[i][1],atoms[i][2],forces[i][0],forces[i][1],forces[i][2])

fordb.write(line)
print "Added "+ocfile+" to force database."

