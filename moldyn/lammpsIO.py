#!/usr/bin/python

from numpy import *
import sys
#mine
from struct_tools import mag,ang,flatten

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
            xlo,xhi,xy=map(float,dump[i+1].split())
            ylo,yhi,xz=map(float,dump[i+2].split())
            zlo,zhi,yz=map(float,dump[i+3].split())
            v1=[xhi-xlo-xy-xz,0,0]
            v2=[xy,yhi-ylo-yz,0]
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


#From a lammps dump read a desired configuration (0-indexed) from the dump file
def dumpReadConfig(dump,configN):
    #Grab the line numbers where the configurations start
    starts=[i for i,line in enumerate(dump) if "TIMESTEP" in line]
    if configN<-1 or configN >= len(starts):
        print "Error: lammpsIO: requested configuration not available in dump file. %d Configurations available."%len(starts)-1
        exit(0)
    if configN==-1:
        configN=len(starts)-1

    #Grab the beginning and ending locations of the configuration in the dump
    start=starts[configN]
    end = len(dump) if len(starts)-1==configN else starts[configN+1]

    #Grab relevant info
    natom=0
    v1,v2,v3=[[0,0,0] for j in range(3)]
    ax,ay,az=[list() for j in range(3)]
    types=list()
    head=""
    i=start
    #File must be parsed in a while loop because there is no garauntee on ordering or that file contains required data
    while i<end:
        line=dump[i]
        if "ITEM:" in line:
            if "TIMESTEP" in line:
                i+=1
                head="From LAMMPS dump, Timestep: %s"%dump[i]
            if "NUMBER OF ATOMS" in line:
                i+=1
                natom=int(dump[i])
                continue

            if "BOX BOUNDS" in line: #assume: ITEM: BOX BOUNDS xy xz yz pp pp p
                i+=1
                xlo,xhi,xy=map(float,dump[i].split())
                i+=1
                ylo,yhi,xz=map(float,dump[i].split())
                i+=1
                zlo,zhi,yz=map(float,dump[i].split())
                v1=[xhi-xlo-xy-xz,0,0]
                v2=[xy,yhi-ylo-yz,0]
                v3=[xz,yz,zhi-zlo]
                continue
            if "ATOMS" in line: #assume ITEM: ATOMS id type x y z
                i+=1
                atominfo=map(lambda x:map(float,x.split()),dump[i:i+natom])
                atominfo.sort(key=lambda x:x[0])
                order,types,ax,ay,az=zip(*atominfo)
                types=map(int,types)
                i+=natom
        else:
            i+=1
    bounds=[v1,v2,v3]
#    if v1[1]!=0 or v1[2]!=0 or v2[0]!=0 or v2[2]!=0 or v3[0]!=0 or v3[1]!=0:
#        print "Warning lammpsIO.readconfig has not been tested for non-orthogonal basis, please double check POSCAR config matches lammps_dump config."
    atoms=zip(ax,ay,az)
    return bounds,types,atoms,head

def dumpWriteConfig(dump,bounds,types,atoms,head):
    data=["#Generated by dumpWriteConfig from dump data\n"]
    data.append( "#"+head.strip("\n").lstrip("#")+"\n\n" )
    data.append( "%d atoms\n"%len(atoms) )
    ntypes=len(set(types))
    data.append( "%d atom types\n\n"%ntypes )
    xhi,yhi,zhi,xy,xz,yz=bounds2lohi(bounds)
    data.append( " 0.0000  % 6.6f  xlo xhi\n"%xhi )
    data.append( " 0.0000  % 6.6f  ylo yhi\n"%yhi )
    data.append( " 0.0000  % 6.6f  zlo zhi\n"%zhi )
    data.append( "% 6.6f  % 6.6f  % 6.6f xy xz yz\n\n"%(xy,xz,yz) )
    data.append( "Atoms\n\n")
    for i,(t,atom) in enumerate(zip(types,atoms)):
        data.append( "\t%d\t%d\t% 6.6f\t% 6.6f\t% 6.6f \n"%(i+1,t,atom[0],atom[1],atom[2]) )
    dump.writelines(data)

#Converts VASP style boundaries to LAMMPS boundaries
def bounds2lohi(bounds):
    v1,v2,v3=map(array,bounds)
    mv1,mv2,mv3=map(mag,[v1,v2,v3]) #magnitude

    origin=array([0,0,0])
    A=ang(v2,origin,v3)
    B=ang(v1,origin,v3)
    C=ang(v1,origin,v2)

    #14 digits of accuracy.
    xhi= mv1
    xy= mv2*cos(C)
    xz= mv3*cos(B)
    yhi=(mv2**2-xy**2)**0.5
    yz= (mv2*mv3*cos(A)-xy*xz)/yhi
    zhi=(mv3**2-xz**2-yz**2)**0.5
    return xhi,yhi,zhi,xy,xz,yz
