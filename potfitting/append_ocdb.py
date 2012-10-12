#!/usr/bin/python
import sys,os
import numpy

#This script writes an OUTCAR configuration to a force database for operating on by potfit

def usage():
    print "Usage:"
    print sys.argv[0].split("/")[-1]+" <force database> <directory containing OUTCAR/POSCAR> <configuration #> <run from script:1>"
    exit(0)

#Arguement locations in sys.argv list
argdb=1
argoc=2
argconfig=3
argscript=4

if not(len(sys.argv) in [4,5]):
    usage()

ocdir = sys.argv[argoc] #outcar directory
outcar = open(ocdir+"/OUTCAR","r") #outcar contains all system information at each timestep/optimation step
poscar = open(ocdir+"/POSCAR","r") #need to use the poscar to fetch the atom types
forcefil= sys.argv[argdb]
fordb = open(forcefil,"r")
grabconfig=int(sys.argv[argconfig])

checkdisable=False
if len(sys.argv)==5:
    if int(sys.argv[argscript])==1:
        checkdisable=True
    else:
        print "If running from a script, make the 4th arguement a \"1\""
        usage()


#get the force configuration number (the total number of configs in the forcefil)
dbcnfgcnt=0 #Force Data Base Configuration Counter
for line in fordb:
    if "#N" in line[0:2]:
        dbcnfgcnt+=1

fordb.close()
fordb = open(forcefil,"a")

poscar.readline()
poscar.readline()
v1=map(float,poscar.readline().split())
v2=map(float,poscar.readline().split())
v3=map(float,poscar.readline().split())
nums=map(int,poscar.readline().split())
types=list()
for i in range(len(nums)):
    types+=[i]*nums[i]
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
    if len(line)==0:
        print "Configuration not found."
        break
    if "FREE ENERGIE OF THE ION-ELECTRON SYSTEM" in line:
        countconfig+=1
        if countconfig==grabconfig:

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

            #Number of atoms, use force, header
            line="#N\t%d 1 ifconf=%d Location:%s  OUTCAR_Config#%d\n"%(natom,dbcnfgcnt,os.getcwd()+"/"+ocdir+"/",countconfig)
            
            #Bounding box (Angstroms)
            line += "#X\t %12.8f  %12.8f %12.8f\n"%(v1[0],v1[1],v1[2])
            line += "#Y\t %12.8f  %12.8f %12.8f\n"%(v2[0],v2[1],v2[2])
            line += "#Z\t %12.8f  %12.8f %12.8f\n"%(v3[0],v3[1],v3[2])

            #Cohesive energy per atom (eV)
            line += "#E\t %12.8f\n"%(PE/natom)

            #Stress Tensor
            stress=[float(s)/1602.0 for s in stresskb.split()] #Convert to floats and divide by 1602kbar to get in proper units
            line += "#S\t %12.8f  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f\n"%(stress[0],stress[1],stress[2],stress[3],stress[4],stress[5]) 
            
            #Positions
            line += "#F\t\n"
            
            for i in range(natom):
                line += "  %d\t %12.8f  %12.8f  %12.8f   %12.8f  %12.8f  %12.8f\n"%(types[i],ax[i],ay[i],az[i],afx[i],afy[i],afz[i])
            #print line
            fordb.write(line)
            print "AWWWW YEAAAAAH. "+ocdir
            break
