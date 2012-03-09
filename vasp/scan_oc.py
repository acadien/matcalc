#!/usr/bin/python
import sys
import numpy

#mine
from duplicate import duplicate26

if len(sys.argv)<2:
    print "Usage:"
    print sys.argv[0]+" <directory containing OUTCAR> <optional: starting config #> <optional: spanning #> <optional 0: disables plotting>"
    exit(0)

outcar = open(sys.argv[1]+"/OUTCAR","r")

poscar = open(sys.argv[1]+"/POSCAR","r") #need to use the poscar to fetch the atom types

poscar.readline()
poscar.readline()
v1=map(float,poscar.readline().split())
v2=map(float,poscar.readline().split())
v3=map(float,poscar.readline().split())
nums=map(int,poscar.readline().split())
types=list()
for i in range(len(nums)):
    types+=[len(nums)-i]*nums[i]

startconfig=0
disabled=False
span=1
if (len(sys.argv)>=3):
    startconfig=int(sys.argv[2])
    if (len(sys.argv)>=4):
        span=int(sys.argv[3])
        if (len(sys.argv)==5):
            if int(sys.argv[4])==0:
                disabled=True

if span==0:
    span=1

if not(disabled):
    from enthought.mayavi import mlab
    fig = mlab.gcf()

N=sum(nums)
count=0
md=True
#NOTE TO SELF: NEVER COMPOSE READ LOOPS LIKE THIS.  TERRIBLE.
while True:
    #Start a new configuration
    line=outcar.readline()
    if len(line)==0:
        break
    if "FORCE on cell" in line:
        count+=1
        if count==startconfig or (count>=startconfig and (count-startconfig)%span==0):
            while True:
                line=outcar.readline()
                if len(line)==0:
                    break
                if "in kB" in line:
                    stresskb=line.split("in kB")[1].strip()
                    break
        
            while True:
                line=outcar.readline()
                if len(line)==0:
                    break
                if "POSITION" in line:
                    outcar.readline()
        
                    ax=list()
                    ay=list()
                    az=list()
                    afx=list()
                    afy=list()
                    afz=list()

                    try:
                        while True:
                            line=outcar.readline()
                            if len(line)==0:
                                break
                            [x,y,z,fx,fy,fz]=map(float,line.split())
                            ax.append(x)
                            ay.append(y)
                            az.append(z)
                            afx.append(fx)
                            afy.append(fy)
                            afz.append(fz)
                    except ValueError:
                        pass
                    break
                
            while True:
                line=outcar.readline()
                if len(line)==0:
                    break
                if "TOTEN" in line:
                    PE=float(line.split("=")[1].split()[0])
                    line=outcar.readline()
                    if md==True:
                        try:
                            KE=float(line.split("=")[1].split()[0])
                            TE=PE+KE
                            T=float(line.split("(temperature")[1].split()[0])
                        except IndexError:
                            md=False
                    break
            print "="*50
            print count
            print "Stress (in kB) XX YY ZZ XY YZ XZ:"
            print stresskb
            if md==True:
                print "T=%g,   PE=%g,   TE=%g" % (T,PE,KE)
            else:
                print "PE=%g" % (PE)
            print "="*50
            print ""
            if not(disabled):
                mlab.clf()
                if md==True:
                    mlab.points3d(ax,ay,az,types,colormap="gist_heat",scale_factor=0.7)
                    mlab.quiver3d(ax,ay,az,afx,afy,afz,line_width=3,scale_factor=1)
                else:
                    datoms,dtypes=duplicate26(zip(ax,ay,az),types,zip(v1,v2,v3))
                    axd,ayd,azd=zip(*datoms)
                    mlab.points3d(axd,ayd,azd,dtypes,colormap="gist_heat",scale_factor=0.7)
                    mlab.plot3d([0,v1[0]],[0,v1[1]],[0,v1[2]],color=(0,1,0))
                    mlab.plot3d([v1[0],v1[0]+v2[0]],[v1[1],v1[1]+v2[1]],[v1[2],v1[2]+v2[2]],color=(0,1,0))
                    mlab.plot3d([v1[0],v1[0]+v3[0]],[v1[1],v1[1]+v3[1]],[v1[2],v1[2]+v3[2]],color=(0,1,0))
                    mlab.plot3d([v1[0]+v2[0],v1[0]+v2[0]+v3[0]],[v1[1]+v2[1],v1[1]+v2[1]+v3[1]],[v1[2]+v2[2],v1[2]+v2[2]+v3[2]],color=(0,1,0))
                
                    mlab.plot3d([0,v2[0]],[0,v2[1]],[0,v2[2]],color=(0,1,0))
                    mlab.plot3d([v2[0],v2[0]+v1[0]],[v2[1],v2[1]+v1[1]],[v2[2],v2[2]+v1[2]],color=(0,1,0))
                    mlab.plot3d([v2[0],v2[0]+v3[0]],[v2[1],v2[1]+v3[1]],[v2[2],v2[2]+v3[2]],color=(0,1,0))
                    mlab.plot3d([v2[0]+v3[0],v1[0]+v2[0]+v3[0]],[v2[1]+v3[1],v1[1]+v2[1]+v3[1]],[v2[2]+v3[2],v1[2]+v2[2]+v3[2]],color=(0,1,0))
                    
                    mlab.plot3d([0,v3[0]],[0,v3[1]],[0,v3[2]],color=(0,1,0))
                    mlab.plot3d([v3[0],v3[0]+v1[0]],[v3[1],v3[1]+v1[1]],[v3[2],v3[2]+v1[2]],color=(0,1,0))
                    mlab.plot3d([v3[0],v3[0]+v2[0]],[v3[1],v3[1]+v2[1]],[v3[2],v3[2]+v2[2]],color=(0,1,0))
                    mlab.plot3d([v3[0]+v1[0],v1[0]+v2[0]+v3[0]],[v3[1]+v1[1],v1[1]+v2[1]+v3[1]],[v3[2]+v1[2],v1[2]+v2[2]+v3[2]],color=(0,1,0))
                mlab.show(stop=True)

print "Total configurations: "+str(count)
