#!/usr/bin/python

import numpy as np
from math import *
from scipy import integrate
#mine
from datatools import superSmooth

#Uses the Ornstein-Zernike relations with the Percus-Yevick closure relationship to
#relate the h(r) back down to the correlation function c(r).

def integ(x,y):
    if len(x)==1:
        return 0
    a=integrate.cumtrapz(y,x=x)
    return a[-1]

def calc_qpr_rltcut(icut,dens,hrx,hry,cr,qr,qrp):
    ti = np.array(range(icut))
    for i,r in enumerate(hrx[:icut]):
        qrp[i] = -r*hry[i] + 2*pi*dens*integ(hrx[:icut],qr[:icut]*(r-hrx[:icut])*hry[np.abs(i-ti)])
    return qrp

def calc_qpr_rgtcut(icut,dens,hrx,hry,cr,qr,qrp):
    N=len(qr)
    qrpLast = np.array(qrp)
    for i,r in enumerate(hrx[icut:]):
        i+=icut
        qrp[i] = -r*cr[i] + 2*pi*dens*integ(hrx[:N-i],qr[:N-i]*qrpLast[i:]) 
    return qrp

def calc_qr(icut,dens,hrx,hry,cr,qr,qrp):
    for i,r in enumerate(hrx):
        qr[i] = -integ(hrx[i:],qrp[i:])
    return qr

def calc_hr_rgtcut(icut,dens,hrx,hry,cr,qr,qrp):
    ti = np.array(range(len(hrx)-icut))+icut
    hryLast = np.array(hry)
    for i,r in enumerate(hrx[icut:]):
        i+=icut
        hry[i] = (-qrp[i] + 2*pi*dens*integ(hrx[icut:],qr[icut:]*(r-hrx[icut:])*hryLast[np.abs(i-ti)]))/r
    return hry

def calc_hr(icut,dens,hrx,hry,cr,qr,qrp):
    ti = np.array(range(len(hrx)))
    hryLast= np.array(hry)
    for i,r in enumerate(hrx):
        if r==0: continue
        hry[i] = (-qrp[i] + 2*pi*dens*integ(hrx,qr*(r-hrx)*hryLast[np.abs(i-ti)]))/r
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

def rdfExtend(rdfX,rdfY,ndens,rmax=50.0,Niter=5,T=100.0,rm=2.5,eps=-1,damped=True,PY=True):
    #rm and eps are parameters for the LJ (used as the PY closure)
    #T is the temperature for the energy, effects smoothing near the transition point
    #PY selects the Percus Yevick Closure for C(r), if false, the HNC criterion are used.
    rdfX=list(rdfX)
    rdfY=list(rdfY)

    drx = rdfX[1]-rdfX[0]
    n = int(np.ceil(rmax/drx))
    hrx = np.array([i*drx for i in range(n)])
    hry = np.array([(j-1.0) for d,j in enumerate(rdfY + [1.0 for i in hrx[len(rdfY):]])])

    phi=[eps*((rm/r)**12-2*(rm/r)**6) if r>0 else eps*((rm/1E-10)**12-2*(rm/1E-10)**6) for r in hrx]
    beta = 1.0/8.6173324e-05/T

    icut = len(rdfX)
    cr = np.array([0.0 for i in hry])
    qr = np.array([0.0 for i in hry])
    qrp = np.array([0.0 for i in hry])

    hrlast = np.array(hry)

    Lmax=10.0
    qbins=1024
    minq,maxq,dq=0,Lmax,Lmax/qbins
    qs=[i*dq+minq for i in range(qbins)]
    qs[0]=1E-10

    import pylab as pl
    for k in range(Niter):
        qrp = calc_qpr_rltcut(icut,ndens,hrx,hry,cr,qr,qrp)
        qrp = calc_qpr_rgtcut(icut,ndens,hrx,hry,cr,qr,qrp)
        qr  = calc_qr        (icut,ndens,hrx,hry,cr,qr,qrp)
        hry = calc_hr_rgtcut (icut,ndens,hrx,hry,cr,qr,qrp)
        cr  = calc_cr_rgtcut (icut,ndens,hrx,hry,cr,qr,qrp,phi,beta,PY)

        err = np.power((hrlast-hry),2).sum()
        print "step %d, error %f"%(k,err)
        if err > 1:
            break
        else:
            temp = np.array(hry)
            calc_hr(len(hrx),ndens,hrx,temp,cr,qr,qrp)
            grExtended = np.array([i+1.0 for i in hry])
            sfExt=list()
            for i,q in enumerate(qs):
                sfExt += [1 + 4*pi*ndens * integrate.simps((grExtended-1.0)*np.sin(q*hrx)*hrx/q ,dx=hrx[1]-hrx[0])]
            #pl.plot(qs,sfExt)
        hrlast = np.array(hry)
    #pl.show()
        
    qrp = superSmooth(hrx,qrp,sigma=0.1)
    calc_hr(len(hrx),ndens,hrx,hry,cr,qr,qrp)
    grExtended = np.array([i+1.0 for i in hry])
    
    return cr,hrx,grExtended
