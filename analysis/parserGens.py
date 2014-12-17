#!/usr/bin/python

#quick and dirty scripts for parsing simple files

def parseEnsemble(ensembleFile,keepAvg=False):
    f = open(ensembleFile,"r")
    f.readline() #eliminate the header
    if keepAvg:
        for line in f:
            try: 
                yield map(float,line.split())
            except ValueError:
                pass
    else:
        for line in f:
            try: 
                yield map(float,line.split()[1:])
            except ValueError:
                pass

def parseOutcarAtoms(byteNums,outcarFile,nAtoms):
    ocar = open(outcarFile,"r")
    
    for i,b in enumerate(byteNums):
        ocar.seek(b)
        ocar.readline()
        ocar.readline()
        atoms = [map(float,ocar.readline().split()[:3]) for a in range(nAtoms)]
        yield atoms

def parseNeighbor(neighborFile):
    f = open(neighborFile,"r")
    f.readline() #eliminate the header
    for line in f:
        try:
            yield map(lambda x: map(int,x.split()),line.split(","))
        except ValueError:
            pass

def parseLammpsBasis(byteNums,lammpsFile):
    f2parse = open(lammpsFile,"r")
    for i,b in enumerate(byteNums):
        f2parse.seek(b)
        f2parse.readline()
        xlo,xhi = map(float,f2parse.readline().split()[:2])
        ylo,yhi = map(float,f2parse.readline().split()[:2])
        zlo,zhi = map(float,f2parse.readline().split()[:2])
        basis = [ [xhi-xlo,0,0], [0,yhi-ylo,0], [0,0,zhi-zlo] ]
        yield basis

def parseLammpsAtoms(byteNums,lammpsFile,nAtoms):
    f2parse = open(lammpsFile,"r")
    for i,b in enumerate(byteNums):
        f2parse.seek(b)
        head = f2parse.readline().split()
        idl = head.index("id")-2
        try:
            ixs = head.index("xs")-2
            iys = head.index("ys")-2
            izs = head.index("zs")-2
        except ValueError:
            ixs = head.index("xu")-2
            iys = head.index("yu")-2
            izs = head.index("zu")-2
        itypes = head.index("type")-2
        atomLines = [f2parse.readline().split() for i in range(nAtoms)]
        atomLines = [map(float,[al[idl],al[ixs],al[iys],al[izs],al[itypes]]) for al in atomLines]
        atomLines.sort()
        ids,ax,ay,az,types  = zip(*atomLines)
        yield zip(ax,ay,az)

def parseLammpsColumns(lammpsFile,nAtoms):
    f2parse = open(lammpsFile,"r")
    headerLen = 9
    seperate = lambda x: map(float,x.split())

    while True:
        lmpData = [f2parse.readline() for i in range(nAtoms + headerLen)]
        yield map(seperate, lmpData[headerLen:])

