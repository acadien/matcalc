#!/usr/bin/python

import sys

def usage():
    print "%s <input force DB> <output force DB> <A=1.0> <B=1.0> <C=0.0>"%sys.argv[0].split("/")[-1]
    print "Re-scales a force database by the ratio of lattice constants"
    print "A = <a_literature / a_adp> = 1.0 -> not altered"
    print "B = sqrt(<Tmelt_literature / Tmelt_adp>) = 1.0 -> not altered"
    print "C = <E_adp - E_literature> = 0.0 -> not altered (ground state energies)"

def readNextEntry(iFDB):

    head = iFDB.pop(0)
    n = int(head.split()[1]) #number of atoms

    bounds = list()
    energy = None
    weights = list()
    stress = list()
    scaleEnable = True
    for i,line in enumerate(iFDB):

        if line[:2] == "#N":
            if "scaleEnable0" in line:
                scaleEnable=False

        if line[:2] == "#C":
            elems = line.split()[1:]

        if line[:2] in ["#X","#Y","#Z"]:
            bounds.append( map( float , line.split()[1:] ) )

        if line[:2] == "#E":
            energy = float(line.split()[1])

        if line[:2] == "#W":
            weights = map( float , line.split()[1:] )

        if line[:2] == "#S":
            stress = map( float , line.split()[1:] )

        if line[:2] == "#F":
            break

    iFDB = iFDB[i+1:]
    atominfo = [ map( float , line.split() ) for line in iFDB[:n] ]
    iFDB = iFDB[n:]

    return [head,elems,bounds,energy,weights,stress,atominfo],iFDB,scaleEnable


def appendEntry(head,elems,bounds,energy,weights,stress,atominfo,oFDB):
    data=[head]
    data.append( "#C\t %s\n"%" ".join(elems) )
    data.append( "#X\t %12.8f  %12.8f %12.8f\n"%tuple(bounds[0]) )
    data.append( "#Y\t %12.8f  %12.8f %12.8f\n"%tuple(bounds[1]) )
    data.append( "#Z\t %12.8f  %12.8f %12.8f\n"%tuple(bounds[2]) )
    data.append( "#E\t "+"%12.8f"%(energy) + "\n" )

    if len(weights)>0:
        data.append( "#W\t "+" ".join(map(str,weights)) + "\n" )

    data.append( "#S\t %12.8f  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f\n"%tuple(stress) )
    data.append( "#F\t\n" )
    data += [ "  %d\t %12.8f  %12.8f  %12.8f   %12.8f  %12.8f  %12.8f\n"%tuple(ai) for ai in atominfo]
    oFDB.writelines(data)

def scaleEntry(fdbEntry,A,B,C):
    head,elems,bounds,energy,weights,stress,atominfo = fdbEntry

    head = head.strip() + " Scale (ABC) %4.4f, %4.4f, %4.4f.\n"%(A,B,C)

    scalePosition = lambda x : x * A
    scaleEnergy = lambda x : x * B - C
    scaleForce = lambda x : x * B / A
    scaleStress = lambda x : x * B / (A*A*A)

    #Scale position, forces and stress by 
    bounds = [ map( scalePosition , i ) for i in bounds] 
    stress = map( scaleStress , stress )
    atominfo = [ [i[0]] + map( scalePosition , i[1:4]) + map( scaleForce , i[4:]) \
                    for i in atominfo ]

    #Scale energy by
    energy = scaleEnergy(energy)

    fdbEntry = [head,elems,bounds,energy,weights,stress,atominfo]
    return fdbEntry

if len(sys.argv) != 6:
    usage()
    exit(0)

iFDB = open(sys.argv[1],"r").readlines()
oFDB = open(sys.argv[2],"w")
aRatio = float(sys.argv[3])
tmRatio = float(sys.argv[4])
eShift = float(sys.argv[5])

entries = sum([1 for line in iFDB if "#N"==line[:2]])

#Loop over all the entries in the force database
for i in range(entries):

    #Grab an entry
    entry, iFDB, scaleEnable = readNextEntry(iFDB)

    #Scale the entry
    if scaleEnable:
        entry = scaleEntry(entry,aRatio,tmRatio,eShift)
    entry.append(oFDB)

    #Stick it in a new force database
    appendEntry(*entry)

