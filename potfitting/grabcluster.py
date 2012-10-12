#!/usr/bin/python
#This script grabs a cluster of atoms from an OUTCAR & POSCAR.  The cluster is within a sphere whose center and radius are defined in the arguements
import sys,thread
from math import sqrt

def plot3(md,ax,ay,az,afx,afy,afz,v1,v2,v3):
    global thelock
    if thelock.locked():
        print "\nPlease close all figures before continuing."
    thelock.acquire()
    mlab.figure()
    if md==True:
        mlab.points3d(ax,ay,az,types,colormap="gist_heat",scale_factor=0.7)
        mlab.axes()
        mlab.quiver3d(ax,ay,az,afx,afy,afz,line_width=3,scale_factor=1)
    else:
        ax2=[0]*N*3
        ay2=[0]*N*3
        az2=[0]*N*3
        ax2[:N]=ax
        ay2[:N]=ay
        az2[:N]=az
        ax2[N:2*N]=[x+v1[0] for x in ax]
        ay2[N:2*N]=[y+v2[1] for y in ay]
        az2[N:2*N]=[z+v3[2] for z in az]
        ax2[2*N:]=[x-v1[0] for x in ax]
        ay2[2*N:]=[y-v2[1] for y in ay]
        az2[2*N:]=[z-v3[2] for z in az]
        mlab.points3d(ax2[:N]*9+ax2[N:2*N]*9+ax2[2*N:]*9,(ay2[:N]*3+ay2[N:2*N]*3+ay2[2*N:]*3)*3,az2*9,types*27,colormap="gist_heat",scale_factor=0.7)
        mlab.axes()
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

    mlab.show()
    thelock.release()
    return


def plot3sphere(md,ax,ay,az,afx,afy,afz,v1,v2,v3,rx,ry,rz,rr):
    global thelock
    if thelock.locked():
        print "\n Please close all other figures before continuing."
    thelock.acquire()

    types2=list()
    for i in types:
        types2.append(i)
    N=len(ax)
    cnt=0
    for i in range(N):
        if sqrt((ax[i]-rx)**2+(ay[i]-ry)**2+(az[i]-rz)**2)<=rr:
            cnt+=1
            if types2[i]==1.0:
                types2[i]=1.25
            else:
                types2[i]=1.75

    print "\nFound %d atoms inside of the sphere." % (cnt)

    if md==True:
        mlab.points3d(ax,ay,az,types2,colormap="gist_heat",scale_factor=0.7)
        mlab.quiver3d(ax,ay,az,afx,afy,afz,line_width=3,scale_factor=1)
    else:
        ax2=[0]*N*3
        ay2=[0]*N*3
        az2=[0]*N*3
        ax2[:N]=ax
        ay2[:N]=ay
        az2[:N]=az
        ax2[N:2*N]=[x+v1[0] for x in ax]
        ay2[N:2*N]=[y+v2[1] for y in ay]
        az2[N:2*N]=[z+v3[2] for z in az]
        ax2[2*N:]=[x-v1[0] for x in ax]
        ay2[2*N:]=[y-v2[1] for y in ay]
        az2[2*N:]=[z-v3[2] for z in az]
        mlab.points3d(ax2[:N]*9+ax2[N:2*N]*9+ax2[2*N:]*9,(ay2[:N]*3+ay2[N:2*N]*3+ay2[2*N:]*3)*3,az2*9,types2*27,colormap="gist_heat",scale_factor=0.7)
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

    sphere = mlab.points3d(rx, ry, rz, scale_mode='none', scale_factor=rr*2, color=(0.67, 0.77, 0.93), resolution=50, opacity=0.7, name='Earth')

    # These parameters, as well as the color, where tweaked through the GUI, with the record mode to produce lines of code usable in a script.
    sphere.actor.property.specular = 0.45
    sphere.actor.property.specular_power = 5
    # Backface culling is necessary for more a beautiful transparent rendering.
    sphere.actor.property.backface_culling = True

    mlab.axes()

    mlab.show()
    thelock.release()
    return

#Reading in the options and preparing for file reading
if len(sys.argv)<2:
    print "Usage:"
    print sys.argv[0]+" <directory containing OUTCAR & POSCAR> <optional: starting config #> <optional 0: disables plotting>"
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
if (len(sys.argv)>=3):
    startconfig=int(sys.argv[2])
    if (len(sys.argv)>=4):
        if int(sys.argv[3])==0:
            disabled=True
if not(disabled):
    from enthought.mayavi import mlab

global thelock
thelock=thread.allocate_lock()

N=sum(nums)
count=0
md=True
started=0
#File reading, plotting and grabbing
while True:
    line=outcar.readline()
    if len(line)==0:
        break
    if "FORCE on cell" in line:
        count+=1
    if count>=startconfig:
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
            print "Step Number: %d"%(count)
            print "Stress (in kB) XX YY ZZ XY YZ XZ:"
            print stresskb
            if md==True:
                print "T=%g,   PE=%g,   TE=%g" % (T,PE,KE)
            else:
                print "PE=%g" % (PE)
            print "="*50
            print ""

            if not(disabled):
                athread=thread.start_new_thread(plot3,(md,ax,ay,az,afx,afy,afz,v1,v2,v3))
            while True:
                try:
                    val=raw_input("Input coords and radius of sphere (x y z r) or (Enter) to continue: ")
                    [rx,ry,rz,rr]=[float(k) for k in val.split()]
                except ValueError:
                    if len(val)==0:
                        break
                else:
                    if not(disabled):
                        bthread=thread.start_new_thread(plot3sphere,(md,ax,ay,az,afx,afy,afz,v1,v2,v3,rx,ry,rz,rr))
                        try:
                            val=raw_input("Good? (y/n):")
                            if val=="y" or val=="Y":
                                N=len(ax)
                                type1=list()
                                type2=list()
                                for i in range(N):
                                    if sqrt((ax[i]-rx)**2+(ay[i]-ry)**2+(az[i]-rz)**2)<=rr:
                                        toprint="% 5.5g  % 5.5g  % 5.5g" % (ax[i],ay[i],az[i])
                                        if types[i]==1:
                                            type1.append(toprint)
                                        else:
                                            type2.append(toprint)
                                print len(type2),len(type1)
                                for i in type2:
                                    print i
                                for i in type1:
                                    print i
                                exit(0)
                        except ValueError:
                            pass
