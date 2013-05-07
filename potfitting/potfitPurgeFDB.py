#!/usr/bin/python
import sys,os
import numpy
#mine
from outcarIO import outcarReadConfig

#This script writes an OUTCAR configuration to a force database for operating on by potfit

def usage():
    print "Usage:"
    print sys.argv[0].split("/")[-1]+" <force database> <configuration #>"
    print "Configuration # is compared against ifconf, if equal, that configuration is removed and all others are decremented."
    exit(0)

if len(sys.argv) != 3:
    usage()

forcefil= sys.argv[1]
fordb = open(forcefil,"r")
grabconfig=int(sys.argv[2])

#get the force configuration number (the total number of configs in the forcefil)
removeFlag = False
nfordb = []
rline = ""
for line in fordb:

    if "#N" in line[0:2]:
        ifconf = int(line.split("ifconf=")[1].split(" ")[0])

        removeFlag = False
        if ifconf==grabconfig:
            rline = str(line)
            removeFlag = True

        if ifconf>grabconfig:
            #Remove 1 from ifconf for configurations after
            pre,post = line.split("ifconf=")
            post = " ".join(post.split()[1:])
            ifconf = ifconf-1
            line = pre + "ifconf=%d "%ifconf + post+"\n"

    if not removeFlag:
        nfordb.append(line)

print ifconf,"configurations"
print "Removed line:"+rline

fordb.close()
open(forcefil,"w").writelines(nfordb)
        
"""
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
print "AWWWW YEAAAAAH. "+ocfile



nums=""
types=list()
ax=list()
ay=list()
az=list()
afx=list()
afy=list()
afz=list()
PE=0
N=sum(nums)
countconfig=0
md=True

#NOTE TO SELF: NEVER COMPOSE READ LOOPS LIKE THIS.  TERRIBLE.
while True:
    #Start a new configuration
    line=outcar.readline()
    if "NIONS" in line and len(nums)==0:
        nums=map(int,line.split("=")[-1].split())        
        for i in range(len(nums)):
            types+=[i]*nums[i]

    if len(line)==0:
        print "Configuration not found."
        break
    if "FREE ENERGIE OF THE ION-ELECTRON SYSTEM" in line:
        countconfig+=1
        if countconfig==grabconfig:
            
            #Lattice vectors
            if "direct lattice vectors" in line:
                line=outcar.readline()
                v1=" ".join(line.split()[0:3])
                line=outcar.readline()
                v2=" ".join(line.split()[0:3])
                line=outcar.readline()
                v3=" ".join(line.split()[0:3])
                break

            #PE
            while True:
                line=outcar.readline()
                if len(line)==0:
                    break
                if "TOTEN" in line:
                    PE=float(line.split("=")[1].split()[0])
                    break

            #STRESS
            while True:
                line=outcar.readline()
                if len(line)==0:
                    break
                if "in kB" in line:
                    stresskb=line.split("in kB")[1].strip()
                    break

            #ATOMIC COORDS & FORCES
            while True:
                line=outcar.readline()
                if len(line)==0:
                    break
                if "POSITION" in line:
                    outcar.readline()
        
                    ax=list()
                    ay=list()
                    az=list()
                    afx=list()
                    afy=list()
                    afz=list()

                    try:
                        while True:
                            line=outcar.readline()
                            if len(line)==0:
                                break
                            [x,y,z,fx,fy,fz]=map(float,line.split())
                            ax.append(x)
                            ay.append(y)
                            az.append(z)
                            afx.append(fx)
                            afy.append(fy)
                            afz.append(fz)
                    except ValueError:
                        pass
                    break
            natom=len(ax)
            
            #Check
            val='n'
            if(not(checkdisable)):
                try:
                    val=raw_input("Write to %s? (yes/no)"%(forcefil))
                except NameError:
                    print "OK not writing to force database, be that way... I'm quitting."
                    break
                if not(val[0] in ["Y","y"]):
                    print "OK not writing to force database, be that way... I'm quitting."
                    break
                
            #
            #Write this configuration to the fordb
            #
"""
