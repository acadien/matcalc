#!/usr/bin/python
from numpy import *
import subprocess
from datatools import flatten

#Returns a POSCAR file contents in string format, just use writeline to write to file or read using poscarIO.read
def outcar2poscar(outcarF,wantConfig=-1):
    basis=zeros([3,3])
    posi=list()

    #Number of atoms from OUTCAR
    grepR = subprocess.check_output("head -n 1000 %s | grep ions\ per\ type "%outcarF,shell=True).split("\n")
    nums = grepR[0].split("=")[1].strip().split()
    Natoms = sum(map(int,nums))

    #Number of timesteps from OUTCAR
    grepR = subprocess.check_output("grep -b Iteration.*1\) %s"%outcarF,shell=True).split("\n")
    grepR = grepR[:-1]
    nsteps = len(grepR)


    #Find the seek location for your desired config
    if wantConfig>nsteps:
        print "Error: outcar2poscar: Requested configuration >%d, the max configurations in OUTCAR."%nsteps
        exit(0)
    seekLocation = int(grepR[wantConfig].split(":")[0])

    wdat="This poscar generated by outcar2poscar.py, from %s configuration #%d\n"%(outcarF,wantConfig)
    wdat+="1.0\n"

    outcar=open(outcarF,"r")
    outcar.seek(seekLocation,0)
    while True:
        #Start a new configuration
        line=outcar.readline()
        if len(line)==0:
            break
        if "FORCE on cell" in line:

            while True:
                line=outcar.readline()
                if "direct lattice vectors" in line:
                    line=outcar.readline()
                    basis[0]=map(float,line.split()[0:3])
                    v1=" ".join(line.split()[0:3])
                    line=outcar.readline()
                    basis[1]=map(float,line.split()[0:3])
                    v2=" ".join(line.split()[0:3])
                    line=outcar.readline()
                    basis[2]=map(float,line.split()[0:3])
                    v3=" ".join(line.split()[0:3])
                    break

            while True:
                line=outcar.readline()
                if "POSITION" in line:
                    outcar.readline()
                    while True:
                        line=outcar.readline()
                        atom=linalg.solve(basis.T,array(map(float,line.split()[0:3])))
                        posi.append(" ".join(map(str,atom)))
                        if len(posi)==Natoms:
                            break
                    break
            break
                
    #Generate the rest of the poscar, return it in the same method as poscarIO.read()
    basis=map(lambda x: map(float,x.split()),[v1,v2,v3])
    atypes = nums
    atoms = map(lambda x:map(float,x.split()),posi)
    head = "Generated by outcarIO.py from %s, configuration %d\n"%(outcarF,wantConfig)
    poscar = ""

    return basis,atypes,atoms,head,poscar


#Returns poscar data and forces in (eV/Angstrom) and stresses (in GPa)
#TE,stress,basis,atoms,forces,types
def outcarReadConfig(outcarF,wantconfig=-1):

    #Ion types and number of steps
    nums = subprocess.check_output("head -n 1000 %s | grep ions\ per\ type "%outcarF,shell=True).split("\n")[0].split("=")[1:]
    nums=map(int,nums)
    types = [[i]*nums[i] for i in range(len(nums))]
    Natoms=sum(nums)

    #Use grep to speed up finding values in a huge file!
    grepResults = subprocess.check_output("grep -b free\ \ energy %s"%outcarF,shell=True).split("\n")
    bytenums=[int(i.split(":")[0]) for i in grepResults if len(i)>2]
    if wantconfig in ["all","All"]:
        wantconfig=range(len(bytenums))

    #Parse the configurations from the specified OUTCAR locations
    if type(wantconfig)==type([]) and len(wantconfig)==1:
        wantconfig=wantconfig[0]
    if type(wantconfig)==int:
        outcar= open(outcarF,"r")
        outcar.seek(bytenums[wantconfig])
        outcar = [outcar.readline() for i in range(5000)]

        TE=float(outcar[0].split("=")[-1].split()[0])

        basis=zeros([3,3])
        atoms=list()
        forces=list()
        stress=list()
        for i,line in enumerate(outcar):
            if "in kB" in line:
                stress=map(lambda x:float(x)/10.0,line.split()[2:])

            if "direct lattice vectors" in line:
                basis[0]=map(float,outcar[i+1].split()[0:3])
                basis[1]=map(float,outcar[i+2].split()[0:3])
                basis[2]=map(float,outcar[i+3].split()[0:3])

            if "POSITION" in line:
                for line in outcar[i+2:i+Natoms+2]:
                    atom=linalg.solve(basis.T,array(map(float,line.split()[0:3])))
                    force=array(map(float,line.split()[3:6]))
                    atoms.append(atom)
                    forces.append(force)
                break
        bt=basis.T
        atoms=array(atoms)
        if sum(atoms[:,0])/len(atoms) < 1.0:
            atomsp=array([bt.dot(atom) for atom in atoms])
        else:
            atomsp=array(atoms)
        return TE,stress,basis,atomsp,forces,types

    elif type(wantconfig)==list:
        TEs=list()
        basiss=list()
        atomss=list()
        forcess=list()
        stresss=list()
        typess=list()
        for wc in wantconfig:
            outcar= open(outcarF,"r") 
            outcar.seek(bytenums[wc])
            oc = [outcar.readline() for i in range(5000)]
            outcar.close()
            outcar=oc

            TE=float(outcar[0].split("=")[-1].split()[0])

            basis=zeros([3,3])
            atoms=list()
            forces=list()
            TE=0
            stress=list()
            for i,line in enumerate(outcar):
                if "in kB" in line:
                    stress=map(lambda x:float(x)/10.0,line.split()[2:])

                if "direct lattice vectors" in line:
                    basis[0]=map(float,outcar[i+1].split()[0:3])
                    basis[1]=map(float,outcar[i+2].split()[0:3])
                    basis[2]=map(float,outcar[i+3].split()[0:3])

                if "POSITION" in line:
                    for line in outcar[i+2:i+Natoms+2]:
                        atom=linalg.solve(basis.T,array(map(float,line.split()[0:3])))
                        force=array(map(float,line.split()[3:6]))
                        atoms.append(atom)
                        forces.append(force)
                    break
            atoms=array(atoms)
            bt=basis.T
            if sum(atoms[:,0])/len(atoms) < 1.0:
                atomsp=array([bt.dot(atom) for atom in atoms])
            else:
                atomsp=array(atoms)

            TEs.append(TE)
            stresss.append(stress)
            basiss.append(basis)
            atomss.append(atomsp)
            forcess.append(forces)
            typess.append(types)
            
        return TEs,stresss,basiss,atomss,forcess,typess

#Returns the number of iterations in an outcar
def nSteps(outcarFile):
    grepR = subprocess.check_output("grep -b Iteration.*1\) %s"%outcarFile,shell=True).split("\n")
    grepR = grepR[:-1]
    n = len(grepR)
    return n

#Returns the size of the timestep in an outcar
def timestep(outcarFile):
    grepR = subprocess.check_output("grep POTIM %s"%outcarFile,shell=True).split("\n")
    ts = float(grepR[0].split()[2])
    return ts

#Returns the number of atoms in the simulation
def nAtoms(outcarFile):
    nums = subprocess.check_output("head -n 1000 %s | grep ions\ per\ type "%outcarFile,shell=True).split("\n")[0].split("=")[1:]
    nums=map(int,nums)
    Natoms=sum(nums)
    return Natoms

#Retursn the types of the atoms in the simulation
def types(outcarFile):
    nums = subprocess.check_output("head -n 1000 %s | grep ions\ per\ type "%outcarFile,shell=True).split("\n")[0].split("=")[1:]
    nums=map(int,nums)
    types = [j for j in flatten([[i]*nums[i] for i in range(len(nums))])]
    return types

#Returns the initial basis of the simulation
def basis(outcarFile):
    spot = int(subprocess.check_output("head -n 5000 %s | grep -b direct\ lattice\ vector"%outcarFile,shell=True).split("\n")[0].split(":")[0])
    
    outcar= open(outcarFile,"r") 
    outcar.seek(spot)
    outcar.readline()
    
    b = [map(float,outcar.readline().split()[0:3]) for i in range(3)]
    
    return b

#Returns the byte wise location of atoms in an OUTCAR file
def atomBytes(outcarFile):
    grepResults = subprocess.check_output("grep -b POSITION %s"%outcarFile,shell=True).split("\n")
    bytenums=[int(i.split(":")[0]) for i in grepResults if len(i)>2]
    return bytenums
