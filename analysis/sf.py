#!/usr/bin/python

from scipy import weave,integrate
from scipy.weave import converters
import numpy as np
from math import *
import cmath
import random
#mine
from rdf import rdf_periodic
from rdfExtend import rdfExtend

#Calculate the Structure factor from the rdf
def sf(rdfX,rdfY,ndens,Lmax=20.0,qbins=1024,damped=True):

    minq,maxq,dq=0,Lmax,Lmax/qbins
    qs=[i*dq+minq for i in range(qbins)]
    qs[0]=1E-10

    maxR = max(rdfX)
    rdfY=np.array(rdfY)
    rdfX=np.array(rdfX)
    dx=rdfX[1]-rdfX[0]

    if damped:
        sf=[1 + 4*pi*ndens * integrate.simps((rdfY-1.0)*np.sin(q*rdfX)*np.sin(pi*rdfX/maxR)*maxR/(q*pi) ,dx=dx) for i,q in enumerate(qs)]
    else:
        sf=[1 + 4*pi*ndens * integrate.simps((rdfY-1.0)*np.sin(q*rdfX)*rdfX/q ,dx=dx) for i,q in enumerate(qs)]

    return qs,sf

def sfq0(rdfX,rdfY,ndens,Lmax=20.0,qbins=1024,damped=None):
    minq,maxq,dq=0,Lmax,Lmax/qbins
    qs=[i*dq+minq for i in range(qbins)]
    qs[0]=1E-10

    maxR = max(rdfX)
    rdfY=np.array(rdfY)
    rdfX=np.array(rdfX)
    dx = rdfX[1]-rdfX[0]

    #Extend the h(r) to get a better estimate near q0
    grExtx,grExty = rdfExtend(rdfX,rdfY,ndens,rmax=50.0,Niter=25,T=1000.0,rm=2.5,eps=-1,damped=True)

    sf1=list()
    for i,q in enumerate(qs):
        sf1 += [1 + 4*pi*ndens * integrate.simps((rdfY-1.0)*np.sin(q*rdfX)*rdfX/q ,dx=dx)]

    import pylab as pl
    
    sf2=list()
    for i,q in enumerate(qs):
        sf2 += [1 + 4*pi*ndens * integrate.simps((grExty-1.0)*np.sin(q*grExtx)*grExtx/q ,dx=dx)]
    #        R=1.0/q
    #        alpha = np.array([(1-rij/2.0/R)**2 * (1+rij/4.0/R) if rij<2*R else 0 for rij in rdfX])
    #        sf2 += [1 + 4*pi*ndens * integrate.simps((rdfY-1.0)*np.sin(q*rdfX)*rdfX/q*alpha,dx=dx)]

    f=open("/home/acadien/Dropbox/sfq0.dat","w")
    a=map(lambda x: "%f %f\n"%(x[0],x[1]),zip(qs,sf2))
    f.writelines(a)
    exit(0)
#    pl.plot(qs,sf1)
#    pl.plot(qs,sf2)
#    pl.show()
#    exit(0)

#stepsize => nqvectors
#0.1 => 425172
#0.2 => 53162
#0.3 => 16138
#0.4 => 7333
#0.5 => 3829
SFbyQCodeGrid="""
double QMAX = qRads[(int)nqbins-1];
double qmag;
double delq = qRads[1]-qRads[0];

double *nqvecs,qdotr,aix,aiy,aiz,qc,qs,dx,dy,dz;
nqvecs = (double*)calloc(nqbins,sizeof(double));
int index;

for(double qx=-QMAX/1.; qx<=QMAX/1.; qx+=qstep){
for(double qy=-QMAX/1.; qy<=QMAX/1.; qy+=qstep){
for(double qz=-QMAX/1.; qz<=QMAX/1.; qz+=qstep){

    qmag=sqrt(qx*qx+qy*qy+qz*qz);
    if(qmag>QMAX) continue;
    if(qmag<0.25) continue;
    index = (int)(qmag/delq);

    qc=0;
    qs=0;
    for(int i=0;i<natoms;i++){
        aix=atoms[i*3+0];
        aiy=atoms[i*3+1];
        aiz=atoms[i*3+2];
        qdotr = qx*aix + qy*aiy + qz*aiz;

        qc += cos(qdotr);
        qs += sin(qdotr);
    }
    sf[index]+= (qc*qc+qs*qs)/natoms;
    nqvecs[index]++;    
}}}

for(int i=0;i<nqbins;i++){
    if(nqvecs[i]<=0.001) continue;
    sf[i]/=nqvecs[i];
}
"""

SFbyQCodeSphere="""
double *qxs,*qys,*qzs,a,b,a2,b2;
int nqs = (int)nqVecs;
qxs = (double*)malloc(sizeof(double)*(int)nqs);
qys = (double*)malloc(sizeof(double)*(int)nqs);
qzs = (double*)malloc(sizeof(double)*(int)nqs);
int qlen=0;
while( qlen < nqs ){
    a = rand()*2/(double)RAND_MAX-1.0;
    b = rand()*2/(double)RAND_MAX-1.0;
    a2 = a*a;
    b2 = b*b;

    if(a2 + b2 > 1.0) 
        continue;

    qxs[qlen] = (2*a*sqrt(1-a2-b2));
    qys[qlen] = (2*b*sqrt(1-a2-b2));
    qzs[qlen] = (1-2*(a2+b2));

    qlen++;
}

double aix,aiy,aiz,qdr;

double r2,*rs;
rs = (double*)malloc(sizeof(double)*natoms);
for(int i=0;i<natoms;i++){
    aix = atoms[i*3];
    aiy = atoms[i*3+1];
    aiz = atoms[i*3+2];
    r2 = aix*aix+aiy*aiy+aiz*aiz;
    rs[i] = sqrt(r2);
}

double qC,qS,qRadius,qx,qy,qz;
double qNorm = 1. / (nqs*natoms);
double cosi,sini;
for(int ql=0;ql<(int)nqbins;ql++){
    qRadius = qRads[ql];

    for(int qv=0;qv<(int)nqs;qv++){
        qx = qxs[qv];
        qy = qys[qv];
        qz = qzs[qv];

        qC=0;
        qS=0;
        for(int i=0;i<natoms;i++){
            aix=atoms[i*3]-basis[0]/2.;
            aiy=atoms[i*3+1]-basis[4]/2.;
            aiz=atoms[i*3+2]-basis[8]/2.;

            qdr = (qx*aix + qy*aiy + qz*aiz)*qRadius;

            qC+=cos(qdr);
            qS+=sin(qdr);
        }
        sf[ql] += (qC*qC+qS*qS);
    }
}

for(int ql=0;ql<(int)nqbins;ql++){
    sf[ql]/=(double)nqs*natoms;
}


"""

def sfq(atoms,basis,nqbins=290,qcut=7.5):

    atoms = np.array([[a[0]%basis[0][0],a[1]%basis[1][1],a[2]%basis[2][2]] for a in atoms])
    natoms = atoms.shape[0]
    basis = np.array(basis)

    #prepare q stuff 
    recip = np.linalg.inv(basis)
    rmodul = map(lambda x: sqrt(sum(x*x)),recip)
    qmin = min(rmodul)
    qRads = np.linspace(qmin,qcut,nqbins)
    sf = np.zeros([nqbins])
    sfsg = np.zeros([nqbins])
    nqVecs = 200
    qstep = 0.1

    #Reshape
    atoms.shape=natoms*3
    basis.shape=9

    #Call the code
    headers=r"""#include <math.h> 
                #include <stdio.h>"""
    weave.inline(SFbyQCodeGrid,['atoms','natoms','qRads','nqbins','nqVecs','sf','basis','sfsg','qstep'],support_code=headers)
    atoms.shape=[natoms,3]
    basis.shape=[3,3]

    import pylab as pl
    pl.plot(qRads,sf,label='tsf')
    pl.plot(qRads,sfsg,label='sfsg')
    pl.legend(loc=0)
    open("blah.SF","w").writelines(["% lf % lf % lf\n"%(q,qc,qs) for q,qc,qs in zip(qRads,sf,sfsg)])
    pl.show()

    return qRads,sf

#number of points = Nz*(Nz-1)+2
def spherePoints(Nz=12):
    #Generate a single circle in the x-y plane
    circleXs = [sin(float(i)/Nz*2*pi) for i in range(Nz)]
    circleYs = [cos(float(i)/Nz*2*pi) for i in range(Nz)]

    #Stagger the circles by some z value
    circleZs = [float(i-Nz/2)/Nz*2.0 for i in range(1,Nz)]

    #Stick them in one list
    xs,ys,zs=[],[],[]
    for z in circleZs:
        r=sqrt(1-z*z)
        xs += [x*r for x in circleXs]
        ys += [y*r for y in circleYs]
        zs += [z for i in range(len(circleXs))]

    #Add top and bottom points to the sphere
    xs = [0] + xs + [0]
    ys = [0] + ys + [0]
    zs = [-1]+ zs + [1]

    return xs,ys,zs

def sphereRandom(Npoints):
    xs,ys,zs=[],[],[]
    while len(xs)<Npoints:
        a,b=random.random()*2-1,random.random()*2-1
        a2=a*a
        b2=b*b
        if a2+b2 >= 1: continue
        
        g=sqrt(1-a2-b2)

        xs.append(2*a*g)
        ys.append(2*b*g)
        zs.append(1-2*(a2+b2))
    """
        xs.append(-2*a*g)
        ys.append(2*b*g)
        zs.append(1-2*(a2+b2))

        xs.append(2*a*g)
        ys.append(-2*b*g)
        zs.append(1-2*(a2+b2))

        xs.append(2*a*g)
        ys.append(2*b*g)
        zs.append(-1+2*(a2+b2))

        xs.append(-2*a*g)
        ys.append(-2*b*g)
        zs.append(1-2*(a2+b2))

        xs.append(-2*a*g)
        ys.append(2*b*g)
        zs.append(-1+2*(a2+b2))

        xs.append(2*a*g)
        ys.append(-2*b*g)
        zs.append(-1+2*(a2+b2))

        xs.append(-2*a*g)
        ys.append(-2*b*g)
        zs.append(-1+2*(a2+b2))
    """     

    return xs,ys,zs

#############################################################################################

SFbyQ2Code="""
double aix,ajx,aiy,ajy,aiz,ajz,c,d;
double rcut2 = rcut * rcut;

//Loop over and pre-calculate the Ri-Rj values.
double *dxs,*dys,*dzs;
int MAXDS = 1000000;
dxs = (double*)malloc(sizeof(double)*MAXDS);
dys = (double*)malloc(sizeof(double)*MAXDS);
dzs = (double*)malloc(sizeof(double)*MAXDS);
int nds=0;

for(int i=0;i<natoms*3;i++){
    aix=atoms[i*3];
    aiy=atoms[i*3+1];
    aiz=atoms[i*3+2];

    for(int j=0;j<natoms*3;j++){
        ajx=atoms[j*3];
        ajy=atoms[j*3+1];
        ajz=atoms[j*3+2];

        //Minimum image distance
        for(int t1=-1;t1<2;t1++){
        for(int t2=-1;t2<2;t2++){
        for(int t3=-1;t3<2;t3++){
            d=aix-ajx+t1*b[0]+t2*b[3]+t3*b[6];
            c=d*d;
            d=aiy-ajy+t1*b[1]+t2*b[4]+t3*b[7];
            c+=d*d;
            d=aiz-ajz+t1*b[2]+t2*b[5]+t3*b[8];
            c+=d*d;

            if(c==0.0) continue;

            if(c<=rcut2){
                dxs[nds] = ajx - aix;
                dys[nds] = ajy - aiy;
                dzs[nds] = ajz - aiz;
                nds++;
            }
        }}}
    }
}

double qC,qS,qRadius;
for(int ql=0;ql<(int)nqbins;ql++){
    qRadius = (double)ql / nqbins * qcut;
    qvals[ql] = qRadius;

    //Loop over q vectors at this radius
    for(int qv=0;qv<(int)nqVecs;qv++){
        qC = 0.0;
        qS = 0.0;
        for(int i=0;i<nds;i++){
            d = (qxs[qv]*dxs[i] + qys[qv]*dys[i] + qzs[qv]*dzs[i])*qRadius;
            qC += cos( d );
            qS += sin( d );
        }
        qbinsCOS[ql]+=qC;
        qbinsSIN[ql]+=qS;
     }
}

double qNorm = 1. / (nqVecs*nds);
for(int ql=0; ql<nqbins;ql++){
    qbinsCOS[ql]*=qNorm;
    qbinsSIN[ql]*=qNorm;
}

free(dxs);
free(dys);
free(dzs);
"""

def sfq2(atoms,basis,nr=100,rcut=12.0,nqbins=100,qcut=4.0):

    atoms = np.array(atoms)
    natoms = atoms.shape[0]
    basis = np.array(basis)

    #Reshape
    atoms.shape=natoms*3
    basis.shape=9
    b=basis.T

    #prepare q stuff 
    qvals = np.zeros([nqbins])
    qbins = np.zeros([nqbins])
    qbinsCOS = np.zeros([nqbins])
    qbinsSIN = np.zeros([nqbins])
    qxs,qys,qzs = map( np.array , spherePoints(1) ) #change the 5 to improve accuracy
    nqVecs = qxs.shape[0]

    #Call the code
    weave.inline(SFbyQ2Code,['atoms','b','natoms','qbinsCOS','qbinsSIN','qvals','nqbins','qcut','rcut','qxs','qys','qzs','nqVecs'])
    atoms.shape=[natoms,3]
    basis.shape=[3,3]

    import pylab as pl
    pl.plot(qvals,qbinsCOS)
    pl.plot(qvals,qbinsSIN)


    #import rdf
    #rdfx,rdfy = rdf.rdf_periodic(atoms,basis)
    #sfx,sfy = sf(rdfx,rdfy,Lmax=6.0)
    #pl.plot(sfx,sfy)
    pl.show()

    return qvals,qbins

#############################################################################################

#weave.inline(SFbyQCode,['atoms','natoms','qvals','qbins','nqbins','qcut','rcut','b','qxs','qys','qzs','nqVecs'])
SFbyQGridCode="""
double aix,ajx,aiy,ajy,aiz,ajz,amx,amy,amz,c,d,di,dj,dx,dy,dz,dij;
double rcut2 = rcut * rcut;
double qC,qS,qRadius;
int nds=0;

for(int i=0;i<natoms;i++){
    aix=atoms[i*3];
    aiy=atoms[i*3+1];
    aiz=atoms[i*3+2];

    for(int j=0;j<natoms;j++){
        if(i==j) continue;

        ajx=atoms[j*3];
        ajy=atoms[j*3+1];
        ajz=atoms[j*3+2];

        dx=ajx-aix;
        dy=ajy-aiy;
        dz=ajz-aiz;
        c=sqrt(dx*dx+dy*dy+dz*dz);

        if(c<=rcut){
            nds++;

            //Loop over q vectors at this radius
            for(int qv=0; qv<nqVecs; qv++){
               di = (qxs[qv]*aix + qys[qv]*aiy + qzs[qv]*aiz);
               //dj = (qxs[qv]*ajx + qys[qv]*ajy + qzs[qv]*ajz);
               //dij = (qxs[qv]*dx + qys[qv]*dy + qzs[qv]*dz);
               qbinsCOS[qrs[qv]] += cos(di);
               qbinsSIN[qrs[qv]] +=1;    
            }
        }
    }
}

//double qNorm = 1. / (nqVecs*natoms);
for(int ql=0; ql<nqbins;ql++){
    qbinsCOS[ql]/=qbinsSIN[ql];
//    qbinsSIN[ql]*=qNorm;
}
"""
