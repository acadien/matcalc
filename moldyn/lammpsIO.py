#!/usr/bin/python

from numpy import *
import sys
import subprocess
import numpy
#mine
from struct_tools import mag,ang

#For reading multiple configurations from a single lammps dump
#Reads the next configuration (jumping by step)
def dumpReadNext(dump,step=0):
    dbreak=0
    chopPoint=0
    for i,line in enumerate(dump):
        if "TIMESTEP" in line:
            dbreak+=1
            if dbreak>step:
                chopPoint=i
                break

    if dbreak==0:
        raise Exception(sys.argv[0].split("/")[-1],'end of file')
    
    dump=dump[chopPoint:]
    end=0
    v1,v2,v3=[[0,0,0] for j in range(3)]
    ax,ay,az=[list() for j in range(3)]
    types=list()
    head=""
    for i,line in enumerate(dump):
        if "TIMESTEP" in line:
            head="From LAMMPS dump, Timestep: %s"%dump[i+1]
            continue
        if "ITEM: BOX BOUNDS" in line:
            if len(dump[i+1].split())==3:
                xlo,xhi,xy=map(float,dump[i+1].split())
                ylo,yhi,xz=map(float,dump[i+2].split())
                zlo,zhi,yz=map(float,dump[i+3].split())
            else:
                xlo,xhi=map(float,dump[i+1].split())
                ylo,yhi=map(float,dump[i+2].split())
                zlo,zhi=map(float,dump[i+3].split())
                xy=xz=yz = 0.0
            v1=[xhi-xlo+xy+xz,0,0]
            v2=[xy,yhi-ylo+yz,0]
            v3=[xz,yz,zhi-zlo]
            continue
        if "ITEM: ATOMS" in line: #assume ITEM: ATOMS id type x y z
            atominfo=list()
            for j,line in enumerate(dump[i+1:]):
                try: 
                    data=map(float,line.split())
                    atominfo.append(data)
                except ValueError: 
                    break
            end=len(atominfo)
            atominfo.sort(key=lambda x:x[0])
            order,types,ax,ay,az=zip(*atominfo)
            types=map(int,types)
            break
    bounds=[v1,v2,v3]
    atoms=zip(ax,ay,az)
    return dump[end:],bounds,types,atoms,head

def parseConfigAtStart(dumpF,seekpoint):
    dumpF.seek(seekpoint)
    while True:
        line = dumpF.readline()
        if "NUMBER OF ATOMS" in line:
            natom=int(dumpF.readline())
            break
    dumpF.seek(seekpoint)

    dump=[dumpF.readline() for i in range(natom+100)]
    v1,v2,v3=[[0,0,0] for j in range(3)]
    ax,ay,az=[list() for j in range(3)]
    types=list()
    head=""
    #File must be parsed in a while loop because there is no garauntee on ordering or that file contains required data
    for i,line in enumerate(dump):
        if "BOX BOUNDS" in line: #assume: ITEM: BOX BOUNDS xy xz yz pp pp p
            if len(dump[i+1].split())==3:
                xlo,xhi,xy=map(float,dump[i+1].split())
                ylo,yhi,xz=map(float,dump[i+2].split())
                zlo,zhi,yz=map(float,dump[i+3].split())
            else:
                xlo,xhi=map(float,dump[i+1].split())
                ylo,yhi=map(float,dump[i+2].split())
                zlo,zhi=map(float,dump[i+3].split())
                xy=xz=yz = 0.0
            v1=[xhi-xlo+xy+xz,0,0]
            v2=[xy,yhi-ylo+yz,0]
            v3=[xz,yz,zhi-zlo]
            continue
        if "ITEM: ATOMS" in line: #assume ITEM: ATOMS id type x y z
            splitfloat = lambda x:map(float,x.split())
            atominfo = map(splitfloat , dump[i+1:i+natom+1])
            atominfo.sort(key=lambda x:x[0])
            order,types,ax,ay,az=zip(*atominfo)[:5]

            if max(ax)<1.1 and min(fabs(ax))>0:
                a,b,c = array(ax),array(ay),array(az)
                ax=v1[0]*a+v2[0]*b+v3[0]*c
                ay=v1[1]*a+v2[1]*b+v3[1]*c
                az=v1[2]*a+v2[2]*b+v3[2]*c

            if v1[1]+v1[2]+v2[0]+v2[2]+v3[0]+v3[1]==0.0:
                delx = v1[0]/2. - sum(ax)/len(ax)
                dely = v2[1]/2. - sum(ay)/len(ay)
                delz = v3[2]/2. - sum(az)/len(az)
                ax = [x+delx for x in ax]
                ay = [y+dely for y in ay]
                az = [z+delz for z in az]

            types=map(int,types)
            break

    basis=[v1,v2,v3]
    atoms=zip(ax,ay,az)

    return basis,types,atoms

#From a lammps dump read a desired configuration (0-indexed) from the dump file
def readConfig(dumpFile,configN):
    #Use grep to speed up finding values in a huge file!
    grepResults = subprocess.check_output("grep -b TIMESTEP %s"%dumpFile,shell=True).split("\n")
    #print dumpFile
    #grepResults = subprocess.check_output("head -n 5 %s"%dumpFile,shell=True).split("\n")
    bytenums=[int(i.split(":")[0]) for i in grepResults if len(i)>2]
    dumpF=open(dumpFile,"r")

    if configN in ["all","All"]:
        configN=range(len(bytenums))

    print "%d configurations in %s"%(len(bytenums),dumpFile)
    try:
        if len(configN)==1:
            configN=configN[0]
    except:
        pass
    if type(configN)==int:
        return parseConfigAtStart(dumpF,bytenums[configN])
    else:
        basiss=list()
        typess=list()
        atomss=list()
        for i in configN:
            b,t,a = parseConfigAtStart(dumpF,bytenums[i])
            basiss.append(b)
            typess.append(t)
            atomss.append(a)
    return basiss,typess,atomss

def dumpWriteConfig(dump,basis,types,atoms,head):
    data=["#Generated by dumpWriteConfig from dump data\n"]
    data.append( "#"+head.strip("\n").lstrip("#")+"\n\n" )
    data.append( "%d atoms\n"%len(atoms) )
    ntypes=len(set(types))
    data.append( "%d atom types\n\n"%ntypes )
    xhi,yhi,zhi,xy,xz,yz=basis2lohi(basis)
    data.append( " 0.0000  % 6.6f  xlo xhi\n"%xhi )
    data.append( " 0.0000  % 6.6f  ylo yhi\n"%yhi )
    data.append( " 0.0000  % 6.6f  zlo zhi\n"%zhi )
    data.append( "% 6.6f  % 6.6f  % 6.6f xy xz yz\n\n"%(xy,xz,yz) )
    data.append( "Atoms\n\n")
    for i,(t,atom) in enumerate(zip(types,atoms)):
        data.append( "\t%d\t%d\t% 6.6f\t% 6.6f\t% 6.6f \n"%(i+1,t,atom[0],atom[1],atom[2]) )
    dump.writelines(data)

#How many atoms in the first step of the dump file
def nAtoms(dump):
    f = open(dump)
    f = [f.readline() for i in range(100)]
    for i,line in enumerate(f):
        if "NUMBER OF ATOMS" in line:
            nAtoms=int(f[i+1])
            break
    return nAtoms

#Returns the location of basis in the dump file
def basisBytes(dump):
    grepResults = subprocess.check_output("grep -b ITEM:\ BOX\ BOUNDS %s"%dump,shell=True).split("\n")
    return [int(i.split(":")[0])+26 for i in grepResults if len(i)>2]

#Returns the location of atoms in the dump file with header
def atomsBytes(dump):
    grepResults = subprocess.check_output("grep -b ITEM:\ ATOMS %s"%dump,shell=True).split("\n")
    return [int(i.split(":")[0]) for i in grepResults if len(i)>2]

#Reads the dump file at the given location and returns a basis set
def parseBasis(dump,b):
    f = open(dump)
    f.seek(b)
    blines=[f.readline() for i in range(3)]
    if len(blines[0].split())==3:
        xlo,xhi,xy=map(float,blines[0].split())
        ylo,yhi,xz=map(float,blines[1].split())
        zlo,zhi,yz=map(float,blines[2].split())
    else:
        xlo,xhi=map(float,blines[0].split())
        ylo,yhi=map(float,blines[1].split())
        zlo,zhi=map(float,blines[2].split())
        xy=xz=yz = 0.0
    v1=[xhi-xlo+xy+xz,0,0]
    v2=[xy,yhi-ylo+yz,0]
    v3=[xz,yz,zhi-zlo]
    return numpy.array([v1,v2,v3])

#Reads the dump file at the given location and returns a list of atoms, nAtoms long
def parseAtoms(dump,b,nAtoms,basis):
    f = open(dump)
    f.seek(b)
    
    #read the header, what index have atom info
    head = f.readline().split()
    ixs = head.index("xs")-2
    iys = head.index("ys")-2
    izs = head.index("zs")-2
    itypes = head.index("type")-2

    #parse the file
    atomLines = [f.readline().split() for i in range(nAtoms)]
    ax,ay,az = zip(*[map(float,[al[ixs],al[iys],al[izs]]) for al in atomLines])
    types = [int(al[itypes]) for al in atomLines]
    
    #move around atoms to fit them into the 
    v1,v2,v3 = basis
    if max(ax)<1.1 and min(fabs(ax))>0:
        a,b,c = array(ax),array(ay),array(az)
        ax=v1[0]*a+v2[0]*b+v3[0]*c
        ay=v1[1]*a+v2[1]*b+v3[1]*c
        az=v1[2]*a+v2[2]*b+v3[2]*c
    if v1[1]+v1[2]+v2[0]+v2[2]+v3[0]+v3[1]==0.0:
        delx = v1[0]/2. - sum(ax)/len(ax)
        dely = v2[1]/2. - sum(ay)/len(ay)
        delz = v3[2]/2. - sum(az)/len(az)
        ax = [x+delx for x in ax]
        ay = [y+dely for y in ay]
        az = [z+delz for z in az]
    
    return numpy.array(zip(ax,ay,az)),numpy.array(types)
    

#Converts VASP style boundaries to LAMMPS boundaries
def basis2lohi(basis):
    v1,v2,v3=map(array,basis)
    mv1,mv2,mv3=map(mag,[v1,v2,v3])
    
    origin=array([0,0,0])
    A=ang(origin,v2,v3) #v2,origin,v3
    B=ang(origin,v1,v3)
    C=ang(origin,v1,v2)

    #14 digits of accuracy.
    xhi=  mv1
    xy =  mv2*cos(C)
    xz =  mv3*cos(B)
    yhi= (mv2**2-xy**2)**0.5
    yz = (mv2*mv3*cos(A)-xy*xz)/yhi
    zhi= (mv3**2-xz**2-yz**2)**0.5

    return xhi,yhi,zhi,xy,xz,yz
