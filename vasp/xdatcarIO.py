#!/usr/bin/python

def read(xdatFile):
    f = open(xdatFile,"r")
    natom = int(f.readline().split()[0])
    [f.readline() for i in range(4)]

    atoms=list()
    count=0
    while True:
        if f.readline() == '':
            break
        atoms.append([map(float,f.readline().split()) for i in range(natom)])
        count+=1

    return atoms
"""
    from mpl_toolkits.mplot3d import Axes3D
    import pylab as pl
    fig = pl.figure()
    ax = fig.gca(projection='3d')
    ax.plot(xs, ys, zs, label='parametric curve')
    pl.show()
    exit(0)
"""
