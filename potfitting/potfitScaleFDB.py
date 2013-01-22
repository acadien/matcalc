#!/usr/bin/python

import sys

def usage():
    print "%s <input force DB> <output force DB> <a_eq/a_adp>"%sys.argv[0].split("/")[-1]
    print "Re-scales a force database by the ratio of lattice constants"
    print "ratio = experimental lattice param / ADP lattice param"

def readNextEntry(iFDB):

    head = iFDB.pop(0)
    n = int(head.split()[1]) #number of atoms
    
    bounds = list()
    energy = None
    weight = None
    stress = list()
    for i,line in enumerate(iFDB):
        
        if line[:2] in ["#X","#Y","#Z"]:
            bounds.append( map( float , line.split()[1:] ) )

        if line[:2] == "#E":
            energy = float(line.split()[1])

        if line[:2] == "#W":
            weight = int(line.split()[1])

        if line[:2] == "#S":
            stress = map( float , line.split()[1:] )

        if line[:2]=="#F":
            break

    iFDB = iFDB[i+1:]
    atominfo = [ map( float , line.split() ) for line in iFDB[:n] ]
    iFDB = iFDB[n:]

    return [head,bounds,energy,weight,stress,atominfo],iFDB


def appendEntry(head,bounds,energy,weight,stress,atominfo,oFDB):
    data=[head]
    data.append( "#X\t %12.8f  %12.8f %12.8f\n"%tuple(bounds[0]) )
    data.append( "#Y\t %12.8f  %12.8f %12.8f\n"%tuple(bounds[1]) )
    data.append( "#Z\t %12.8f  %12.8f %12.8f\n"%tuple(bounds[2]) )
    data.append( "#E\t "+"%12.8f"%(energy) + "\n" )

    if weight!=None:
        data.append( "#W\t "+str(int(weight)) + "\n" )

    data.append( "#S\t %12.8f  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f\n"%tuple(stress) )
    data.append( "#F\t\n" )
    data += [ "  %d\t %12.8f  %12.8f  %12.8f   %12.8f  %12.8f  %12.8f\n"%tuple(ai) for ai in atominfo]
    oFDB.writelines(data)


def scaleEntry(fdbEntry,ratio):
    head,bounds,energy,weight,stress,atominfo = fdbEntry

    head = head.strip() + " Scaled by %5.5f.\n"%ratio

    scaleMul = lambda x : x*ratio
    scaleDiv = lambda x : x/ratio

    #Scale position and force data by 
    bounds = [ map( scaleMul , i ) for i in bounds] 
    stress = map( scaleDiv , stress )
    atominfo = [[i[0]] + map( scaleMul , i[1:]) for i in atominfo]

    fdbEntry = [head,bounds,energy,weight,stress,atominfo]
    return fdbEntry

if len(sys.argv) != 4:
    usage()
    exit(0)

iFDB = open(sys.argv[1],"r").readlines()
oFDB = open(sys.argv[2],"w")
ratio= float(sys.argv[3])

entries = sum([1 for line in iFDB if "#N"==line[:2]])
for i in range(entries):
    entry, iFDB = readNextEntry(iFDB)

    entry = scaleEntry(entry,ratio)

    entry.append(oFDB)
    appendEntry(*entry)
