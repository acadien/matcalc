#!/usr/bin/python

import sys
import operator
from numpy import *
import scipy
#mine
from poscarIO import readposcar

#Multiplies by volume to get absolute Charge
def readchgcar(chgcar):
    v1,v2,v3,types,xs,ys,zs,header,chgcar=readposcar(chgcar)
    poscardata=(v1,v2,v3,types,xs,ys,zs,header)

    chgcar.pop(0)
    gridsz=[int(i) for i in chgcar.pop(0).split()]

    Tot_pnts = reduce(operator.mul,gridsz)
    vol=dot(v1,cross(v2,v3))/Tot_pnts

    #*vol
    chgdata=array([float(i) for i in "".join(chgcar).split()[:Tot_pnts]])
    chgdata=asarray(chgdata)
    chgdata.shape=gridsz
    chgdata=swapaxes(chgdata,0,2)

    return poscardata,gridsz,chgdata

#Writes a chargecar to vtk format
def writeVTK(fname,poscardata,gridsz,chgdata,clean=False):
    v1,v2,v3,types,xs,ys,zs,header=poscardata
    
    Tot_pnts=reduce(operator.mul,gridsz)

    header=list()
    header.append("# vtk DataFile Version 2.0")
    header.append("Charge Density 3D grid")
    header.append("ASCII")
    header.append("DATASET STRUCTURED_POINTS")
    header.append("DIMENSIONS "+" ".join(map(str,gridsz)))
    header.append("ORIGIN 0.0 0.0 0.0")
    header.append("SPACING "+" ".join(map(str,[i/j for i,j in zip([v1[0],v2[1],v3[2]],gridsz)])))
    header.append("POINT_DATA %d"%Tot_pnts)
    header.append("SCALARS charge_density float")
    header.append("LOOKUP_TABLE default")
    
    vtkfile=open(fname,"w")
    vtkfile.writelines([x+"\n" for x in header])

    if clean==True:
        chgdata=chgdata.ravel()
        def clearoutlier(i,a): return 1.0 if i>a else i/a
        avgval=scipy.sum(chgdata)/Tot_pnts
        vclearout=vectorize(clearoutlier)
        chgdata=vclearout(chgdata,avgval)
        #chgdata=array([i/avgval for i in chgdata])
        chgdata.shape=gridsz
    #chgdata=chgdata.ravel()
    for i in range(gridsz[0]):
         for j in range(gridsz[1]):
            vtkfile.write(" ".join(map(lambda x:str(x),chgdata[i][j]))+"\n")


def usage():
    print "Usage: %s <input chgcar> <optional:output vtk filename>"%sys.argv[0].split("/")[-1]

if len(sys.argv)<2:
    usage()
    exit(0)
    
chgfname=sys.argv[1]
if len(sys.argv)==3:
    vtkfname=sys.argv[2]
else:
    vtkfname=chgfname+".vtk"

poscardata,gridsz,chgdata=readchgcar(open(chgfname,"r").readlines())
writeVTK(vtkfname,poscardata,gridsz,chgdata,clean=True)