#!/usr/bin/python


#From a lammps dump read a desired configuration (0-indexed) from the dump file
def dumpReadConfig(dump,configN):
    #Grab the line numbers where the configurations start
    starts=[i for i,line in enumerate(dump) if "TIMESTEP" in line]
    if configN < 0 or configN > len(starts):
        print "Error: lammpsIO: requested configuration not available in dump file. %d Configurations available."%len(starts)-1
        exit(0)

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
                v1=[xhi-xlo,0,0]
                v2=[xy,yhi-ylo,0]
                v3=[xz,yz,zhi-zlo]
                continue
            if "ATOMS" in line: #assume ITEM: ATOMS id type x y z
                i+=1
                types,ax,ay,az=zip(*map(lambda x:map(float,x.split()[1:]),dump[i:i+natom]))
                types=map(int,types)
                i+=natom
        else:
            i+=1
    bounds=[v1,v2,v3]
    if v1[1]!=0 or v1[2]!=0 or v2[0]!=0 or v2[2]!=0 or v3[0]!=0 or v3[1]!=0:
        print "Warning lammpsIO.readconfig has not been tested for non-orthogonal basis, please double check POSCAR config matches lammps_dump config."
    atoms=zip(ax,ay,az)
    return bounds,types,atoms,head
