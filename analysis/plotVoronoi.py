#!/usr/bin/python
import sys
from numpy import *
import pylab as pl
import mpl_toolkits.mplot3d.art3d as art3d
import mpl_toolkits.mplot3d.axes3d as p3
#mine
from voronoiIO import *
from poscarIO import readposcar
from poscarsPlot import plotsimulation
from struct_tools import *

def usage():
    print "Usage: %s <qvoronoi o-file (Qo_output)> <qvoronoi Fi-file (QFi_output)> <POSCAR>"%(sys.argv[0])
    
def plot_polyhedron(polys,ax=None):
    #Polys is a multidimensional list such that:
    #polys[i]=a list/array of polygons corresponding to each polyhedron
    #polys[i][j]=a list/array of 3-D points corresponding to each polygon
    if ax==None:
        fig=pl.figure()
        ax = p3.Axes3D(fig)

    poly3d=art3d.Poly3DCollection(polys)#,facecolor="none",edgecolor="black")   
    ax.add_collection3d(poly3d)

if __name__=="__main__":
    if len(sys.argv)<4:
        usage()
        exit(0)

    qvodata=open(sys.argv[1],"r").readlines()
    qvfidata=open(sys.argv[2],"r").readlines()
    poscar=open(sys.argv[3],"r").readlines()

    fig=pl.figure()
    ax3d = p3.Axes3D(fig)

    [v1,v2,v3,atypes,ax,ay,az,head,poscar] = readposcar(poscar)
    basis=array([v1,v2,v3])
    atoms=zip(ax,ay,az)
    bounds=array([basis[0][0],basis[1][1],basis[2][2]])

    polyverts=readQVo(qvodata)#,bounds)
    polyplanes,neighbors=readQVFi(qvfidata)

    #Chop off polyhedra of points outside simulation
    polyplanes=polyplanes[:len(atoms)]

    polyhedra=[points2polyhedron(verts,planes,plotting=True) for verts,planes in zip(polyverts,polyplanes)]
    #polyhedra=zip(polyverts,polyplanes)
    #plotsimulation(basis,atoms,atypes,ax3d)
        
    #plot_polyhedron(polyhedra,ax3d)

    ax3d.set_xlabel('X')
    ax3d.set_ylabel('Y')
    ax3d.set_zlabel('Z')

    print len(polyhedra)

    for poly in polyhedra[:108]:
        plot_polyhedron(poly,ax3d)

    ax3d.set_xlim3d([-5,20])
    ax3d.set_ylim3d([-5,20])
    ax3d.set_zlim3d([-5,20])
    pl.show()

