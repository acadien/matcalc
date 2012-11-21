#!/usr/bin/python
#This script grabs a cluster of atoms from an OUTCAR & POSCAR.  The cluster is within a sphere whose center and radius are defined in the arguements
import sys,threading
from math import sqrt
from os import path
#mine
from poscarIO import readposcar
from datatools import flatten

def plot3(md,ax,ay,az,afx,afy,afz,v1,v2,v3,types):
    N=len(ax)
    global thelock
    if thelock.isSet():
        print "\nPlease close all figures before continuing."
    thelock.set()
    mlab.figure()
    if md==True:
        mlab.points3d(ax,ay,az,[i+1 for i in types],colormap="gist_heat",scale_factor=0.1,vmin=-1.0,vmax=4.0)#colormap="gist_heat"
        mlab.axes()
        mlab.quiver3d(ax,ay,az,afx,afy,afz,line_width=3,scale_factor=1)
    else:
        print "here"
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
        print len(types)*27,len(az2)*9
        mlab.points3d(ax2[:N]*9+ax2[N:2*N]*9+ax2[2*N:]*9,(ay2[:N]*3+ay2[N:2*N]*3+ay2[2*N:]*3)*3,az2*9,[d+1.0 for d in types]*27,scale_factor=1.0,vmin=-1.0,vmax=4.0)
        #mlab.points3d(ax,ay,az,[d+1.0 for d in types],scale_factor=1.0)
        mlab.axes()
    mlab.plot3d([0,v1[0]],[0,v1[1]],[0,v1[2]],color=(0,1,0),line_width=0.1)
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

#    mlab.show(stop=True)
#    engine=mlab.get_engine()
    mlab.show(stop=True)
    thelock.clear()
    return


def plot3sphere(md,ax,ay,az,afx,afy,afz,v1,v2,v3,rx,ry,rz,rr):
    global thelock
    if thelock.isSet():
        print "\n Please close all other figures before continuing."
    thelock.set()

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
    mlab.clf()
    if md==True:
        mlab.points3d(ax,ay,az,types,scale_factor=1,vmin=-1.0,vmax=4.0)
#        mlab.quiver3d(ax,ay,az,afx,afy,afz,line_width=3,scale_factor=1)
    else:
        mlab.points3d(ax,ay,az,types,scale_factor=1,vmin=-1.0,vmax=4.0)#colormap="gist_heat"
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
        mlab.points3d(ax2[:N]*9+ax2[N:2*N]*9+ax2[2*N:]*9,(ay2[:N]*3+ay2[N:2*N]*3+ay2[2*N:]*3)*3,az2*9,types2*27,scale_factor=1.0,vmin=-1.0,vmax=4.0)
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

    sphere = mlab.points3d(rx, ry, rz, scale_mode='none', scale_factor=rr*2, color=(0.67, 0.77, 0.93), resolution=50, opacity=1.0, name='Earth')

    # These parameters, as well as the color, where tweaked through the GUI, with the record mode to produce lines of code usable in a script.
    sphere.actor.property.specular = 0.45
    sphere.actor.property.specular_power = 5
    # Backface culling is necessary for more a beautiful transparent rendering.
    sphere.actor.property.backface_culling = True

    mlab.axes()

    mlab.show(stop=True)
    thelock.clear()
    return

def prompt_iface(md,count,stresskb,T,PE,KE,ax,ay,az,afx,afy,afz,v1,v2,v3,types):
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

    if not(plotdisabled):
        athread=threading.Thread(target=plot3,args=(md,ax,ay,az,afx,afy,afz,v1,v2,v3,types))
        athread.start()
        athread.join()
    while True:
        try:
            val=raw_input("Input coords and radius of sphere (x y z r) or (Enter) to continue:\n")
            [rx,ry,rz,rr]=[float(k) for k in val.split()]
        except ValueError:
            if len(val)==0:
                break
        else:
            if not(plotdisabled):
                bthread=threading.Thread(target=plot3sphere,args=(md,ax,ay,az,afx,afy,afz,v1,v2,v3,rx,ry,rz,rr))
                bthread.start()
                bthread.join()
                try:
                    val=raw_input("Good? (y/n):")
                    if val=="y" or val=="Y":

                        N=len(ax)
                        atoms=[[ax[i],ay[i],az[i],types[i]] for i in range(N) if sqrt((ax[i]-rx)**2+(ay[i]-ry)**2+(az[i]-rz)**2)<=rr]

                        atoms.sort(key=lambda x:x[3])
                        d1,d2,d3,types=zip(*atoms)
                        cs=[types.count(i) for i in range(N)]
                        #cut off types list after last non-zero number
                        last=0
                        for i,v in enumerate(cs):
                            if v>0:
                                last=i
                        cs=cs[:last+1]
                        print cs[:]
                        for atom in atoms:
                            print "% 5.5g  % 5.5g  % 5.5g" % (tuple(atom[:3]))
                            
                        exit(0)
                except ValueError:
                    pass

def usage():
    print "Usage:"
    print "%s <POSCAR> <XDATCAR> <optional: starting config #> <optional 0: disables plotting>"%(sys.argv[0].split("/")[-1])

#Reading in the options and preparing for file reading
if len(sys.argv)<3:
    usage()
    exit(0)

    
poscar = open(sys.argv[1],"r")
xdatcar = open(sys.argv[2],"r").readlines()[6:]

v1,v2,v3,atypes,ax,ay,az,head,poscar=readposcar(poscar.readlines())
afx,afy,afz=[0]*len(ax),[0]*len(ay),[0]*len(az)
types=[i for i in flatten([[i]*num for i,num in enumerate(atypes)])]
startconfig=0
plotdisabled=False

if (len(sys.argv)>=4):
    startconfig=int(sys.argv[3])
if (len(sys.argv)>=5):
    if int(sys.argv[4])==0:
        plotdisabled=True
if not(plotdisabled):
    from mayavi import mlab

global thelock
thelock=threading.Event()

count=0
md=True
started=0
stresskb=T=PE=KE=0

#File reading, plotting and grabbing
Natoms=len(ax)
for i in range(startconfig,int(len(xdatcar)/Natoms)):
    ax,ay,az=zip(*[map(float,line.split()) for line in xdatcar[i*(1+Natoms):(i+1)*(Natoms+1)-1]])
    
    prompt_iface(md,count,stresskb,T,PE,KE,ax,ay,az,afx,afy,afz,[1,0,0],[0,1,0],[0,0,1],types)
    """
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
            prompt_iface(md,count,stresskb,T,PE,KE,ax,ay,az,afx,afy,afz,v1,v2,v3,types)
    """
