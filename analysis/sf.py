#!/usr/bin/python

from scipy import weave
from scipy.weave import converters
import scipy
import numpy as np
from math import *
import cmath

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
        sf=[1 + 4*pi*ndens * scipy.integrate.simps((rdfY-1.0)*np.sin(q*rdfX)*np.sin(pi*rdfX/maxR)*maxR/(q*pi) ,dx=dx) for i,q in enumerate(qs)]
    else:
        sf=[1 + 4*pi*ndens * scipy.integrate.simps((rdfY-1.0)*np.sin(q*rdfX)*rdfX/q ,dx=dx) for i,q in enumerate(qs)]

    return qs,sf

#weave.inline(SFbyQCode,['atoms','natoms','qvals','qbins','nqbins','qcut','rcut','b','qxs','qys','qzs','nqVecs'])
SFbyQCode="""
double aix,ajx,aiy,ajy,aiz,ajz,amx,amy,amz,c,d,di,dj,dx,dy,dz,dij;
double rcut2 = rcut * rcut;
double qC,qS,qRadius;
int nds=0;

for(int ql=0;ql<(int)nqbins;ql++){
    qRadius = (double)ql / nqbins * qcut;
    qvals[ql] = qRadius;
}

for(int i=0;i<natoms*3;i++){
    aix=atoms[i*3];
    aiy=atoms[i*3+1];
    aiz=atoms[i*3+2];

    for(int j=0;j<natoms*3;j++){
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
/*
            for(int ql=0;ql<(int)nqbins;ql++){
                qRadius = qvals[ql];
                if(qRadius==0.0) qRadius=1E-10;
                if(c==0) err=1;
                qbinsCOS[ql]+=sin(qRadius*c)/qRadius/c;
            }
*/
            for(int ql=0;ql<(int)nqbins;ql++){
                qRadius = qvals[ql];

                if(qRadius == 0.0)
                   qRadius+=1E-5;

                //Loop over q vectors at this radius
                for(int qv=0; qv<nqVecs; qv++){
                    //di = (qxs[qv]*aix + qys[qv]*aiy + qzs[qv]*aiz)*qRadius;
                    //dj = (qxs[qv]*ajx + qys[qv]*ajy + qzs[qv]*ajz)*qRadius;
                    dij = (qxs[qv]*dx + qys[qv]*dy + qzs[qv]*dz)*qRadius;
                    qbinsCOS[ql] += cos(dij);
                    //qbinsSIN[ql] += (cos( di ) * cos( dj ) + sin( di ) * sin( dj ));
                }
            }
        }
    }
}

//double qNorm = 1. / (nqVecs*natoms);
//for(int ql=0; ql<nqbins;ql++){
//    qbinsCOS[ql]*=qNorm;
//    qbinsSIN[ql]*=qNorm;
//}
"""

def sfq(atoms,basis,nr=100,rcut=10.0,nqbins=200,qcut=10.0):

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
    qxs,qys,qzs = map( np.array , spherePoints(5) ) #change the 5 to improve accuracy
    nqVecs = qxs.shape[0]

    err = 0
    #Call the code
    weave.inline(SFbyQCode,['atoms','b','natoms','qbinsCOS','qbinsSIN','qvals','nqbins','qcut','rcut','qxs','qys','qzs','nqVecs','err'])
    atoms.shape=[natoms,3]
    basis.shape=[3,3]

    print max(qbinsCOS),max(qbinsSIN)
    print err
    qbinsCOS/=max(qbinsCOS)
    #qbinsSIN/=max(qbinsSIN)
    import pylab as pl
    pl.plot(qvals,qbinsCOS)
    #pl.plot(qvals,qbinsSIN)



    #import rdf
    #rdfx,rdfy = rdf.rdf_periodic(atoms,basis)
    #sfx,sfy = sf(rdfx,rdfy,Lmax=6.0)
    #pl.plot(sfx,sfy)
    pl.show()

    return qvals,qbins

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

def gridPoints(Nz=4):
    gp=list()
    Nz1=Nz+1
    Nz=float(Nz)
    for i in range(Nz1):
        for j in range(Nz1):
            for k in range(Nz1):
                gp.append([i/Nz,j/Nz,k/Nz])
    return zip(*gp)

#number of points = Nz*(Nz-1)+2
def spherePoints(Nz=12):
    #Generate a single circle in the x-y plane
    circleXs = [sin(float(i)/Nz*2*pi) for i in range(Nz)]
    circleYs = [cos(float(i)/Nz*2*pi) for i in range(Nz)]
    
    #Stagger the circles by some z value
    circleZs=[float(i-Nz/2)/Nz*2.0 for i in range(1,Nz)]

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

    xs,ys,zs = zip(*[[x,y,z] for x,y,z in zip(xs,ys,zs) if x>=0 and y>=0 and z>=0])

    return xs,ys,zs

