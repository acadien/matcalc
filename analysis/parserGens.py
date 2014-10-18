#!/usr/bin/python

#quick and dirty scripts for parsing simple files

def parseEnsemble(ensembleFile):
    f = open(ensembleFile,"r")
    f.readline() #eliminate the header
    for line in f:
        try: 
            yield map(float,line.split()[1:])
        except ValueError:
            pass

def parseOutcarAtoms(bytenums,outcarFile,nAtoms):
    for i,b in enumerate(bytenums):
        outcarFile.seek(b)
        outcarFile.readline()
        outcarFile.readline()
        atoms = [map(float,outcarFile.readline().split()[:3]) for a in range(nAtoms)]
        yield atoms

def parseNeighbor(neighborFile):
    f = open(neighborFile,"r")
    f.readline() #eliminate the header
    for line in f:
        try:
            yield map(lambda x: map(int,x.split()),line.split(","))
        except ValueError:
            pass
