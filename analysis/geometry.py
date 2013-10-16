#!/usr/bin/python

from scipy import *
#mine
from struct_tools import *

#Given a list of atoms and periodic boundary conditions, returns a list of pairs of indeces that are atoms at that distance
#key part: fabs(dist(a,b)-length) < err => true
def atomsAtLength(atoms,halfNeighbors,length,err=0.05,periodic=False,bounds=None):
    keypairs=list()
    for iNeighb,Neighbs in enumerate(halfNeighbors):
        atomi=atoms[iNeighb]
        for jNeighb in Neighbs:
            atomj=atoms[jNeighb]
            if periodic:
                if fabs(dist_periodic(atomi,atomj,bounds)-length) < err:
                    keypairs.append([iNeighb,jNeighb])
            else:
                if fabs(dist(atomi,atomj)-length) < err:
                    keypairs.append([iNeighb,jNeighb])

    return keypairs

def cutPeriodicBounds(pnt,bnds):
    code = """
    int tot=0;
    for(int i=0;i<3;i++){
        if(pnt[i]>=bnds[i]){
            pnt[i]=bnds[i];
            tot++;
        }
        if(pnt[i]<0){
            pnt[i]=0;
            tot++;
        }
    }
    if(tot==3)
        return_val=-1;
    else
        return_val=0;
    """
    if weave.inline(code,['pnt','bnds']) == -1:
        return None
    return pnt

def applyPeriodicBounds(pnt,bnds):
    code = """
    for(int i=0;i<3;i++){
        while(pnt[i]>=bnds[i])
            pnt[i]-=bnds[i];
        while(pnt[i]<0)
            pnt[i]+=bnds[i];
    }
    """
    weave.inline(code,['pnt','bnds'])      
    return pnt

#Return a set of points on the plane defined by the miller indeces.
def miller2pnts(mInd,N=5):
    #mind = [a,b,c] miller indeces
    #N = number of points on line, N^2 points are returned.
    
    [1./m for m in mInd if m>0]
        

#===========  Find Structures ==================
def find_squares(atoms,neighbs,atol=1.0,dtol=0.1):#tolerance in degrees
    atol=radians(atol)
    N=len(atoms)
    right=pi/2.0
    squares=list()
    uniqs=list()

    for i,ai in enumerate(atoms):
        for j in neighbs[i]:
            aj=atoms[j]
            dij=dist(ai,aj)
            for k in neighbs[i]:
                if j==k: continue
                ak=atoms[k]
                angle = ang(ai,aj,ak)
                if fabs(dist(aj,ak)-dij) > dtol: continue
                if fabs(angle - right) <= atol:
                    fourth = set(neighbs[k]) - set([i,j,k])
                    for l in fourth:
                        al = atoms[l]
                        if fabs(dist(ai,al)-dij) > dtol: continue
                        if fabs(dist(ak,al)-dij) > dtol: continue
                        square=[i,j,k,l]
                        ssq = sorted(square)
                        if ssq not in uniqs:
                            squares.append(square)
                            uniqs.append(ssq)

    return squares
     
def find_cubes(squares):
    multiverts=[v for v in flatten(squares)]
    vertices=[v for v in set(multiverts)]

    #These vertices are a part of 3 different squares 
    rightverts = [v for v in vertices if multiverts.count(v) >= 3]

    #Check if each value is in some biglist
    def allin(values,biglist):
        for value in values:
            if value not in biglist:
                return False
        return True
    
#    for i,isq in enumerate(squares):
#        if allin(isq,rightverts):
#            print isq


#================= Making Polyhedron from points (Voronoi) ============

#Given a list of points all in the same hyperplane, defined by the coefs from 'hyperplane' return the list of points sorted and ready to be used as a polygon.
def points2polygon(points,hyperplane):
    com=array([sum(i)/len(i) for i in zip(*points)])
    points2=[i-com for i in points]
    s=min([i for i in zip(range(len(points2)),points2) if i[1][2]>0],key=lambda x:x[1][2])
    rot=rotmatx(points2[0],array([1,0,0]))
    points2=[dot(rot,p) for p in points2]
    com=zeros(3)
    start=points2[0]
    angles=[ang(start,com,i) for i in points2]
    for i,p in enumerate(points2):
        if p[2]<0 and points2[0][2]>0:
            angles[i]=2*pi-angles[i]
        if p[2]>0 and points2[0][2]<0:
            angles[i]=2*pi-angles[i]
    p=zip(angles,points)
    p.sort(key=lambda x:x[0])
    return zip(*p)[1]


#Given a list of points that belong to the same polyhedron and
#the list of hyperplanes that those points are on
#Return a list of polygons that make that polyhedron
def points2polyhedron(points,planes,plotting=True):
    polyhedron=[list() for i in range(len(planes))]
    for i,plane in enumerate(planes):
        cnt=0
        for j,point in enumerate(points):
            a=fabs(dot(point,plane[0:3])+plane[3])
            if a<1e-9:
                cnt+=1
                polyhedron[i].append(point)
        if cnt==0:
            print "ERROR: points2polyhedron, hyperplane with no vertices... something is wrong here."
            exit(0)
        #Need to sort out the polygons if you plan on plotting these
        if plotting:
            polyhedron[i]=points2polygon(polyhedron[i],plane)
    return polyhedron
