#!/usr/bin/python

from scipy import *
import operator
#mine
import poscarIO


#poscardata,chgdata = chgcarIO.read(open(chgcarFileName,"r").readlines())
#if spin polarized, returns poscardata,chgdata_add,chgdata_sub
#sum(sum(sum(chgdata))) = Total # of electrons (or valence electrons)
def read(chgcar,SP=False,frac_coord=True):
    basis,types,atoms,header,chgcar=poscarIO.read(chgcar,frac_coord)
    poscardata=(basis,types,atoms,header)

    chgcar.pop(0)
    gridsz=[int(i) for i in chgcar.pop(0).split()]

    Tot_pnts = reduce(operator.mul,gridsz)

    npline=len(chgcar[0].split())
    cur=Tot_pnts/npline
    if cur*npline<Tot_pnts:
        cur+=1
    #Divices by Tot_pnts to get absolute Charge (e.g. total # electrons)
    chgdata=array(map(lambda x:float(x)/Tot_pnts,"".join(chgcar[:cur]).split()))
    chgdata.shape=gridsz
    chgcar=chgcar[cur:]

    #Spin polarized chgcar, first half is chgdata-summed, second is subtracted
    if SP:
        if len(chgcar)>0:
            chgdata_add=chgdata
            chgcar = chgcar[2:]
            chgdata_sub=array(map(lambda x:float(x)/Tot_pnts,"".join(chgcar[:cur]).split()))
            chgdata_sub.shape=gridsz
        else:
            print "WARNING: Requested Spin-Polarized CHGCAR read from a non-spin polarized CHGCAR file."
            SP=False

    #To eliminate points that step outside a reasonable range
    avg = sum(chgdata)/len(chgdata)

    #chgdata[where(chgdata > 3*avg)[0]]=avg


    #Some chgcars seem to need their axes swapped, possible bug in VASP
#    chgdata=swapaxes(chgdata,0,2)
    if SP:
        return poscardata,chgdata_add,chgdata_sub
    return poscardata,chgdata

#Writes a chargecar to vtk format
def writeVTK(fname,poscardata,chgdata):
    basis,types,atoms,header=poscardata
    gridsz=chgdata.shape

    Tot_pnts=reduce(operator.mul,gridsz)

    header=list()
    header.append("# vtk DataFile Version 2.0")
    header.append("Charge Density 3D grid")
    header.append("ASCII")
    header.append("DATASET STRUCTURED_POINTS")
    header.append("DIMENSIONS "+" ".join(map(str,gridsz)))
    header.append("ORIGIN 0.0 0.0 0.0")
    header.append("SPACING "+" ".join(map(str,[i/j for i,j in zip([basis[0][0],basis[1][1],basis[2][2]],gridsz)])))
    header.append("POINT_DATA %d"%Tot_pnts)
    header.append("SCALARS charge_density float")
    header.append("LOOKUP_TABLE default")
    
    vtkfile=open(fname,"w")
    vtkfile.writelines([x+"\n" for x in header])

    mx=chgdata.max()
#    mx=1
    for i in range(gridsz[0]):
         for j in range(gridsz[1]):
            vtkfile.write(" ".join(map(lambda x:str(x/mx),chgdata[i][j]))+"\n")
