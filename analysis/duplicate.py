#!/usr/bin/python

from scipy import *

#duplicates the atoms in the 6 (+/-) axial, 12 planar and 8 translations of those bases, making a 3x3x3 cube
def duplicate26(atoms,types,basis):
    v1,v2,v3=zip(*basis)
    [ax,ay,az]=zip(*atoms)
    N=len(types)

    ax=list(ax)
    ay=list(ay)
    az=list(az)
    
    #Shift atoms to center ("lowest" atom is by 0,0,0)
    #ax=[x+v1[0]+v2[0]+v3[0] for x in ax]
    #ay=[y+v1[1]+v2[1]+v3[1] for y in ay]
    #az=[z+v1[2]+v2[2]+v3[2] for z in az]

    #Make new atomic coordinates
    axp=[x+v1[0]+v2[0]+v3[0] for x in ax]
    ayp=[y+v1[1]+v2[1]+v3[1] for y in ay]
    azp=[z+v1[2]+v2[2]+v3[2] for z in az]
    axm=[x-v1[0]-v2[0]-v3[0] for x in ax]
    aym=[y-v1[1]-v2[1]-v3[1] for y in ay]
    azm=[z-v1[2]-v2[2]-v3[2] for z in az]

    #Make new basis
    v1=[i*3 for i in v1]
    v2=[i*3 for i in v2]
    v3=[i*3 for i in v3]

    #Duplicate
    dax = ax*9+axp*9+axm*9
    day = ay*3+ayp*3+aym*3
    day *= 3
    daz = az+azp+azm
    daz *= 9

    datoms = zip(dax,day,daz)
    dtypes = types*27
    dbasis = zip(v1,v2,v3)

    return array(datoms),dtypes,dbasis



#duplicates the atoms along the 3 (+) basis directions, and 4 diagonals, making a 2x2x2 cube
def duplicate7(atoms,types,basis):
    v1,v2,v3=zip(*basis)
    [ax,ay,az]=zip(*atoms)
    N=len(types)

    ax=list(ax)
    ay=list(ay)
    az=list(az)

    #Make new atomic coordinates
    axp=[x+v1[0]+v2[0]+v3[0] for x in ax]
    ayp=[y+v1[1]+v2[1]+v3[1] for y in ay]
    azp=[z+v1[2]+v2[2]+v3[2] for z in az]

    #Make new basis
    v1=[i*2 for i in v1]
    v2=[i*2 for i in v2]
    v3=[i*2 for i in v3]
    
    #Duplicate
    dax = ax*4+axp*4
    day = ay*2+ayp*2
    day *= 2
    daz = az+azp
    daz *= 4

    datoms = zip(dax,day,daz)
    dtypes = types*8
    dbasis = zip(v1,v2,v3)

    return datoms,dtypes,dbasis
