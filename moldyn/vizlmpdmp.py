#!/usr/bin/python
import sys
from enthought.mayavi import mlab
import pylab as pl

def min2(a,b):
    if abs(a)<abs(b):
        return a
    return b

rcut=2.644877

if len(sys.argv)==1:
    print "Usage: ./plotatoms.py <positions file> [optional: skip to] [optional:step-span]"# [optional:dump]"
    print "The 3rd, optional, arguement is the number at which to start the viewing."
    print "The 3rd, optional, arguement is the number of positions to skip when plotting the atomic configuration.  Negative values are ignored."
    #print "If the 4th argument is given (in addition to the 3rd), any value is fine, than the option to dump one of the configurations to the configuration database is given."
    exit(0)

span=0
skip=0
if len(sys.argv)>=3:
    if int(sys.argv[2])>=0:
        skip=int(sys.argv[2])
if len(sys.argv)>=4:
    if int(sys.argv[3])>=0:
        span=int(sys.argv[3])

xs=list()
ys=list()
zs=list()
ts=list()
times=[0]
bndx=[0]*3
bndy=[0]*3
bndz=[0]*3
Natoms=0

cnt=0
with open(sys.argv[1],"r") as afil:
    #Grab the header info just once.
    while True:
        if "ITEM: NUMBER OF ATOMS" in afil.readline():
            Natoms = int(afil.readline())
            break

    while True:
        if "ITEM: BOX BOUNDS" in afil.readline():
            bndx = map(float,afil.readline().split())
            bndy = map(float,afil.readline().split())
            bndz = map(float,afil.readline().split())
            lx=bndx[1]
            ly=bndy[1]
            lz=bndz[1]
            break

    #Grab the atom locations for as many times its in there...
    while True:
        line = afil.readline()
        if len(line)==0: 
            break
        if "ITEM: TIMESTEP" in line:
            times.append(float(afil.readline()))
        if "ITEM: ATOMS" in line:
            xs.append(list())
            ys.append(list())
            zs.append(list())
            ts.append(list())
            for i in range(Natoms):
                [a,t,x,y,z]=map(float,afil.readline().split())
                xs[-1].append(x)
                ys[-1].append(y)
                zs[-1].append(z)
                ts[-1].append(t)
        cnt+=1

#Calculate the pair correlation function.
Npnts=1000
rcut=max(xs[0])
width=30
xdummy=[[float(i)*rcut for i in range(Npnts)]]*width
ydummy=[[j]*Npnts for j in range(width)]
pairdistro=list()

for x,y,z in zip(xs,ys,zs):
    radii=list()
    pairdistro.append([[0]*Npnts]*width)
    for i in range(Natoms):
        for j in range(i+1,Natoms):
            radii.append(((x[i]-x[j])**2+(y[i]-y[j])**2+(z[i]-z[j])**2)**0.5)
    
    for r in radii:
        if r<rcut:
            pairdistro[-1][0][int(r*Npnts/rcut)]+=1

    norm=float(sum(pairdistro[-1][0]))/100000.0
    for i,a in enumerate(pairdistro[-1][0]):
        a=float(a)/norm
        for j in range(width):
            pairdistro[-1][j][i]=a

#print pairdistro[0][0],sum(pairdistro[0][0])
    
fig1 = mlab.figure()
fig2 = mlab.figure()
#mlab.gcf().scene.y_plus_view()
view2=mlab.view(120,90)
print view2
#mlab.pitch(90)
span+=1
for j in range(skip,cnt-1):
    if j%span!=0 and j!=0:
        continue
    print "TimeStep: %d" % (times[j])
    xx=xs[j]
    yy=ys[j]
    zz=zs[j]
    tt=ts[j]
    mlab.clf(fig1)
    mlab.clf(fig2)
    mlab.points3d(xx,yy,zz,tt,colormap="gist_heat",scale_factor=0.2,figure=fig1)
    mlab.mesh(xdummy,ydummy,pairdistro[j],figure=fig2)
    mlab.show(stop=True)

mlab.close(all=True)


#Bond lengths
"""
print "Calculating histogram of bond lengths corresponding to the minimum energy."
n=pes.index(min(pes))
xs=xs[n][:N]
ys=ys[n][:N]
zs=zs[n][:N]
radii=list()
for i in range(N):
    for j in range(i+1,N):
        radii.append(((xs[i]-xs[j])**2+(ys[i]-ys[j])**2+(zs[i]-zs[j])**2)**0.5)
bins=pl.array(range(1000))/1000.0*rcut*1.5
pl.hist(radii,bins)
pl.show()
"""
"""
#Writing the optimal configuration to a file.
dbdir="./configDB/"
if(len(sys.argv)==4):
    print "Which position would you like to store in the configDB?[optimal]"
    try:
        n=int(input())
        fname=sys.argv[1].split(".")[0]+".cfg"
        print "Writing the "+str(n)+" position & energy to the configuration file "+fname
    except SyntaxError:
        n=pes.index(min(pes))
        fname=sys.argv[1].split(".")[0]+".cfg"
        print "Writing the (optimal)"+str(n)+" position & energy to the configuration file "+fname
    
    cfgfil=open(dbdir+fname,"w")
    cfgfil.write("#Configuration file for "+str(N)+" atoms using the Dzugutov potential"+"\n")
    
    cfgfil.write("POTENTIAL ENERGY : "+str(pes[n])+"\n")
    cfgfil.write("NUM ATOMS        : "+str(N)+"\n")
    cfgfil.write("BOND ENERGY      : "+str(pes[n]/N)+"\n")
    cfgfil.write("POSITIONS        : "+"\n")
    for i in range(N):
        dt="% 7.7e % 7.7e % 7.7e\n" % (xs[i],ys[i],zs[i])
        cfgfil.write(dt)
    #maybe do some calculations on the positions (bond lengths, pair distribution)
    cfgfil.write("HISTOGRAM OF BOND LENS : \n")
    radii=list()
    for i in range(N):
        for j in range(i+1,N):
            radii.append(((xs[i]-xs[j])**2+(ys[i]-ys[j])**2+(zs[i]-zs[j])**2)**0.5)
    bins=pl.array(range(1000))/1000.0*rcut*1.5
    binned,dummy,dummy2=pl.hist(radii,bins)
    cfgfil.write(str(len(binned))+"\n")
    for i in range(len(binned)):
        dt="% 7.7e %d\n" % (bins[i]+rcut*1.5/2000,binned[i])
        cfgfil.write(dt)
"""
