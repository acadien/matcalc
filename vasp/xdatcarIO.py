#!/usr/bin/python
import subprocess

def read(xdatFile):
    startBlockFlag=False
    headlen=0
    nAtomsTot = 0
    atoms = [[]]
    for i,line in enumerate(open(xdatFile,"r")):
        if i==0:
            nAtomsTot = int(line.split()[0])
            continue

        if len(line.split())==0:
            headlen=i
            startBlockFlag=True
            continue

        if startBlockFlag:
            try:
                atoms[-1].append(map(float,line.split()))
            except ValueError:
                atoms.pop(-1)
                atoms.append([])
                startBlockFlag=False
            if len(atoms[-1])==nAtomsTot:
                atoms.append([])
                startBlockFlag=False

    atoms.pop(-1)
    print "%d configs found in %s"%(len(atoms),xdatFile)
    return atoms

def read2(xdatFile):
    startBlockFlag=False
    headlen=0
    nAtomsTot = 0
    atoms = []

    f=open(xdatFile,"r")
    nAtomsTot = int(f.readline().split()[0])
    [f.readline() for i in range(4)]
    atomRange = range(nAtomsTot)

    while f.readline():
        atoms = map(lambda x: map(float,x.split()),f.read(nAtomsTot*38).split("\n")[:-1])
        #atoms = [map(float,f.readline().split()) for a in atomRange]
        yield atoms

def nSteps(xdatFile):
    natom = int(open(xdatFile,"r").readline().split()[0])
    grepR = subprocess.check_output("wc -l %s"%xdatFile,shell=True).split("\n")
    n = int(grepR[0].split()[0])/(natom+1)
    return n
"""
    from mpl_toolkits.mplot3d import Axes3D
    import pylab as pl
    fig = pl.figure()
    ax = fig.gca(projection='3d')
    ax.plot(xs, ys, zs, label='parametric curve')
    pl.show()
    exit(0)
"""
