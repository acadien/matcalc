#!/usr/bin/python

from math import *
import numpy as np
#mine
from theta2q import d2theta,theta2q,d2q,q2d

# alpha=gamma=90

#HKLs selected for alpha'' structure for Cerium
monohkls=["20-1","110","200","11-1","002","20-2","111","201","31-1","11-3","31-2","202","020","310","40-2","31-3","20-4","311","220","022","113","22-2","004"]

def prephkl(hkl):
    seperate=list()
    i=0
    while i<len(hkl):
        if hkl[i]=="-":
            seperate.append(float(hkl[i:i+2]))
            i+=1
        else:
            seperate.append(float(hkl[i]))
        i+=1
    return seperate

def calcd(a,b,c,beta,hkl):
    h,k,l=prephkl(hkl)
    ssb=sin(radians(beta))**2
    cb=cos(radians(beta))
    return 1./sqrt( (h**2/a**2 + ssb*k**2/b**2 + l**2/c**2 - 2*h*l*cb/a/c)/ssb )

def calcds(a,b,c,beta,hkls=None):
    if hkls==None:
        hkls=monohkls
    return [calcd(a,b,c,beta,hkl) for hkl in hkls]

def calc2theta(a,b,c,beta,lamda,hkl):
    return d2theta(calcd(a,b,c,beta,hkl),lamda)

def calc2thetas(a,b,c,beta,lamda,hkls=None):
    return map(lambda g:d2theta(g,lamda),calcds(a,b,c,beta,hkls))

def calcq(a,b,c,beta,hkl):
    return d2q(calcd(a,b,c,beta,hkl))

def calcqs(a,b,c,beta,hkls=None):
    return map(d2q,calcds(a,b,c,beta,hkls))

def struct_lsq(qs,hkls=None):
    if hkls==None:
        hkls=monohkls
    hs,ks,ls=map(np.array,zip(*map(prephkl,hkls)))
    A=np.array([hs**2,ks**2,ls**2,-2*hs*ls]).T
    y=np.array(map(q2d,qs))**(-2)
    (coeffs, residuals, rank, sing_vals) = np.linalg.lstsq(A, y)
    v0,v1,v2,v3 = coeffs
    beta=acos(v3*(v0*v2)**(-0.5))
    a=(v0*sin(beta)**2)**(-0.5)
    c=(v2*sin(beta)**2)**(-0.5)
    b=v1**(-0.5)
    beta=degrees(beta)
    print "(a,b,c,beta) = ( %4.4g, %4.4g, %4.4g, %4.4g )"%(a,b,c,beta)
    print "residual = ",residuals
    print map(lambda x: round(x,3),calcqs(a,b,c,beta))
    print monohkls
