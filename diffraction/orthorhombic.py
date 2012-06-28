#!/usr/bin/python

from math import *
import numpy as np
#mine
from theta2q import d2theta,theta2q,d2q,q2d

# alpha=beta=gamma=90
orthohkls=["020","110","002","021","111","022","112","130","131","200","023","040","113","041","132","220","202","221","004","042"]

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

def calcd(a,b,c,hkl):
    h,k,l=prephkl(hkl)
    return 1.0/sqrt( h**2  / a**2 + k**2 / b**2 + l**2 / c**2)

def calcds(a,b,c,hkls=None):
    if hkls==None:
        hkls=orthohkls
    return [calcd(a,b,c,hkl) for hkl in hkls]

def calc2theta(a,b,c,lamda,hkl):
    return d2theta(calcd(a,b,c,hkl),lamda)

def calc2thetas(a,b,c,lamda,hkls=None):
    return map(lambda g:d2theta(g,lamda),calcds(a,b,c,hkls))

def calcq(a,b,c,hkl):
    return d2q(calcd(a,b,c,hkl))

def calcqs(a,b,c,hkls=None):
    return map(d2q,calcds(a,b,c,hkls))

def struct_lsq(qs,hkls=None):
    if hkls==None:
        hkls=orthohkls
    hs,ks,ls=map(np.array,zip(*map(prephkl,hkls)))
    A=np.array([hs**2,ks**2,ls**2]).T
    y=np.array(map(q2d,qs))**(-2)
    (coeffs, residuals, rank, sing_vals) = np.linalg.lstsq(A, y)
    print "(a,b,c) = ( %4.4g, %4.4g %4.4g )"%tuple(np.array(coeffs)**(-0.5))
    print "residual = ",residuals
    print map(lambda x: round(x,4),calcqs(*(coeffs**(-0.5))))
    print orthohkls
