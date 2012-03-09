#!/usr/bin/python
import sys
from numpy import matrix,linalg

def poscarFixScale(thefile,prompt=1):
    fil=open(thefile,"r")
    points=list()
    cnt=0
    lc=0
    out=""
    vecs=list()
    for line in fil:
        lc+=1
        if lc<6 and lc>=3:
            vecs.append([float(i) for i in line.split()])
        if lc<8:
            out+=line
            continue
        if len(line)<=1: continue #skip newlines
        [x,y,z]=[float(i) for i in line.split()]
        if x<=1 and y<=1 and z<=1:
            cnt+=1
        points.append([x,y,z])

    if cnt==len(points):
        print "This POSCAR appears to already be scaled.  Exiting."
        exit(0)

    fil.close()

    A = matrix(vecs)

    for p in points:
        x=linalg.solve(A,p)[:]
        out += "% 5.9f % 5.9f % 5.9f\n" % (x[0],x[1],x[2])

    if prompt==1:
        try:
            raw_input("Overwriting file \"%s\" with newly rescaled positions. Hit (enter) to continue."%thefile)
        except NameError:
            pass


    open(thefile,"w").write(out)
    
if __name__=="__main__":
    if len(sys.argv)!=2:
        print "Usage:"
        print "%s <poscar directory> "%sys.argv[0]
        exit(0)
    fixposcarscale(sys.argv[1].strip("/")+"/POSCAR")
