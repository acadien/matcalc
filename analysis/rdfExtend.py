#!/usr/bin/python

import numpy as np
from math import *
from scipy import integrate
#mine
from datatools import superSmooth
from colors import vizSpec

#Uses the Ornstein-Zernike relations with the Percus-Yevick closure relationship to
#relate the h(r) back down to the correlation function c(r).

def integ(x,y):
    if len(x)==1:
        return 0
    a=integrate.cumtrapz(y,x=x)
    return a[-1]

def calc_qpr_rltcut(icut,dens,hrx,hry,cr,qr,qrp):
    ti = np.array(range(len(hrx)))
    for i,r in enumerate(hrx[:icut]):
        qrp[i] = -r*hry[i] + 2*pi*dens*integ(hrx,   qr * (r-hrx) * hry[np.abs(i-ti)]    )
    return qrp

def calc_qpr_rgtcut(icut,dens,hrx,hry,cr,qr,qrp):
    N=len(qr)
    qrpLast = np.array(qrp)
    for i,r in enumerate(hrx[icut:]):
        i+=icut
        qrp[i] = -r*cr[i] + 2*pi*dens*integ( hrx[:N-i] ,qr[:N-i]*qrpLast[i:]) 
        #qrp[i] = -r*cr[i] + 2*pi*dens*integ( hrx[:N-i] ,qr[:N-i]*qr[i:]) 
        
    return qrp

def calc_qr(icut,dens,hrx,hry,cr,qr,qrp):
    for i,r in enumerate(hrx):
        qr[i] = -integ(hrx[i:], qrp[i:])
    return qr

def calc_hr_rgtcut(icut,dens,hrx,hry,cr,qr,qrp):
    ti = np.array(range(len(hrx)))
    hryLast = np.array(hry)
    for i,r in enumerate(hrx[icut:]):
        i+=icut
        hry[i] = (-qrp[i] + 2*pi*dens*integ(hrx,qr * (r-hrx) * hryLast[np.abs(i-ti)]) ) /r
    return hry

def calc_cr_rgtcut(icut,dens,hrx,hry,cr,qr,qrp,phi,beta,PY):
    gr =hry+1.0
    for i,r in enumerate(hrx[icut:]):
        i+=icut
        if PY: #Use percus-yevick
            cr[i] = (1.0-np.exp(phi[i]*beta))*gr[i]
        else: #Use HyperNettedChain HNC
            cr[i] = hry[i]-np.log(gr[i]*np.exp(phi[i]*beta))
    return cr

def calc_cr(icut,dens,hrx,hry,cr,qr,qrp,phi,beta,PY):
    gr =hry+1.0
    for i,r in enumerate(hrx):
        if r>1:
            if PY:
                cr[i] = (1.0-np.exp(phi[i]*beta))*gr[i]
            else:
                cr[i] = hry[i]-np.log(gr[i]*np.exp(phi[i]*beta))
    return cr

def calcSq(ndens,hrx,cr,qs):
    dx = hrx[1]-hrx[0]
    sfcq=list()
    for i,q in enumerate(qs):
        sfcq.append( 1.0/(1.0 - 4.*pi*ndens/q * integrate.simps(cr*np.sin(q*hrx)*hrx, dx=dx)))
    return sfcq

def calcSqh(ndens,hrx,hry,qs):
    dx = hrx[1]-hrx[0]
    sfExt=list()
    for i,q in enumerate(qs):
        sfExt += [1 + 4*pi*ndens/q * integrate.simps(hry*np.sin(q*hrx)*hrx ,dx=dx)]
    return sfExt

def calc_hr(icut,dens,hrx,hry,cr,qr,qrp):
    ti = np.array(range(len(hrx)))
    hryLast= np.array(hry)
    hry2 = np.array(hry)
    for i,r in enumerate(hrx):
        if r==0: continue
        hry2[i] = (-qrp[i] + 2*pi*dens*integ(hrx,qr*(r-hrx)*hryLast[np.abs(i-ti)]))/r
    return hry2

"""
def calc_hr(icut,dens,hrx,hry,cr,qr,qrp):
    ti = np.array(range(len(hrx)))
    hryLast= np.array(hry)
    for i,r in enumerate(hrx):
        if r==0: continue
        hry[i] = (-qrp[i] + 2*pi*dens*integ(hrx,qr*(r-hrx)*hryLast[np.abs(i-ti)]))/r
    return hry

def calc_cr(icut,dens,hrx,hry,cr,qr,qrp,phi,beta,PY):
    gr =hry+1.0
    for i,r in enumerate(hrx):
        if r>1:
            if PY:
                cr[i] = (1.0-np.exp(phi[i]*beta))*gr[i]
            else:
                cr[i] = hry[i]-np.log(gr[i]*np.exp(phi[i]*beta))
    return cr
"""

def rdfExtend(rdfX,rdfY,ndens,rmax=50.0,Niter=5,T=300.0,rm=2.5,epsilon=-1E-10,damped=True,PY=True):
    #rm and eps are parameters for the LJ (used as the PY closure)
    #T is the temperature for the energy, effects smoothing near the transition point
    #PY selects the Percus Yevick Closure for C(r), if false, the HNC criterion are used.
    rdfX=list(rdfX)[:-10]
    rdfY=list(rdfY)[:-10]

    drx = rdfX[1]-rdfX[0]
    n = int(np.ceil(rmax/drx))
    hrx = np.array([i*drx for i in range(n)])
    hry = np.array([(j-1.0) for d,j in enumerate(rdfY + [1.0 for i in hrx[len(rdfY):]])])

    phi=[ epsilon*((rm/r)**12-2*(rm/r)**6) if r>0 else 0.0 for r in hrx] #else eps*((rm/1E-10)**12-2*(rm/1E-10)**6) for r in hrx]
    #import pylab as pl
    #pl.plot(hrx,phi)
    #pl.ylim(-10,100)
    #pl.show()
    #exit(0)
    beta = 1.0/8.6173324e-05/T

    icut = len(rdfX)
    cr = np.array([0.0 for i in hry])
    qr = np.array([0.0 for i in hry])
    qrp = np.array([0.0 for i in hry])

    hrlast = np.array(hry)
    qrplast = np.zeros(len(hry))
    Lmax=10.0
    qbins=1024
    minq,maxq,dq=0,Lmax,Lmax/qbins
    qs=[i*dq+minq for i in range(qbins)]
    qs[0]=1E-10

    PY=True
    eta = 0.2

    qs=[i/100.0 for i in range(1,1000)]

    import pylab as pl
    kbreak=20
    for k in range(Niter):
        if k==kbreak: break
        qrp = calc_qpr_rltcut(icut,ndens,hrx,hry,cr,qr,qrp)
        qrp = calc_qpr_rgtcut(icut,ndens,hrx,hry,cr,qr,qrp)
        qrp = eta*qrp + (1-eta)*qrplast
        qr  = calc_qr        (icut,ndens,hrx,hry,cr,qr,qrp)
        hry = calc_hr_rgtcut (icut,ndens,hrx,hry,cr,qr,qrp)
        cr  = calc_cr (icut,ndens,hrx,hry,cr,qr,qrp,phi,beta,PY)
        
        hry2 = calc_hr (icut,ndens,hrx,hry,cr,qr,qrp)
        sq = calcSqh(ndens,hrx,hry2,qs)
        #pl.plot(hrx,hry,label=str(k),c=vizSpec(float(k)/kbreak))
        pl.plot(qs,sq,label=str(k),c=vizSpec(float(k)/kbreak))
        err = np.power((hrlast-hry),2).sum()
        print "step %d, error %f"%(k,err)
        if err>100: break

        hrlast = np.array(hry)
        qrplast = np.array(qrp)


    pl.legend(loc=0)
    pl.show()
    exit(0)
    qrp = superSmooth(hrx,qrp,sigma=0.1)
    calc_hr(len(hrx),ndens,hrx,hry,cr,qr,qrp)
    pl.plot(hrx,hry,c="black",lw=2)
    grExtended = np.array([i+1.0 for i in hry])

    return cr,hrx,grExtended
