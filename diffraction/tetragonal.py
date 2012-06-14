#!/usr/bin/python

from math import *
import numpy as np
#mine
from theta2q import d2theta,theta2q,d2q,q2d

# a=b, alpha=beta=gamma=90

def calcd(a,c,hkl):
    h,k,l=map(float,zip(*hkl)[0])
    return 1.0/sqrt((h**2 + k**2) / a**2 + l**2 / c**2)

def calcds(a,c,hkls=None):
    if hkls==None:
        hkls=["101","002","110","112","200","103","211","202","004"]
    return [calcd(a,c,hkl) for hkl in hkls]

def calc2theta(a,c,lamda,hkl):
    return d2theta(calcd(a,c,hkl),lamda)

def calc2thetas(a,c,lamda,hkls=None):
    return map(lambda g:d2theta(g,lamda),calcds(a,c,hkls))

def calcq(a,c,hkl):
    return d2q(calcd(a,c,hkl))

def calcqs(a,c,hkls=None):
    return map(d2q,calcds(a,c,hkls))

def struct_002_110(q002,q110):
    d110=q2d(q110)
    d002=q2d(q002)
    c=sqrt(4.0*d002**2)
    a=sqrt(2.0*d110**2)
    print "(a,c) = ( %4.4g, %4.4g )"%(a,c)
    print map(lambda x: round(x,4),calcqs(a,c))

def struct_lsq(qs,hkls=None):
    if hkls==None:
        hkls=["101","002","110","112","200","103","211","202","004"]
    hs,ks,ls=[np.array(map(float,g)) for g in zip(*hkls)]
    A=np.array([hs**2+ks**2,ls**2]).T
    y=np.array(map(q2d,qs))**(-2)
    (coeffs, residuals, rank, sing_vals) = np.linalg.lstsq(A, y)
    print "(a,c) = ( %4.4g, %4.4g )"%tuple(np.array(coeffs)**(-0.5))
    print "residual = ",residuals
    print map(lambda x: round(x,4),calcqs(*(coeffs**(-0.5))))
