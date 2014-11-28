#!/usr/bin/python

from math import *

#Returns Qs
#def theta2q(2theta,lamda):
#    return [4.0*pi*sin(i/2.0)/lamda for i in 2theta] 
def theta2q(theta,lamda):
    return 4.0*pi*sin(radians(theta/2.0))/lamda

def q2theta(q,lamda):
    return 2.0*degrees(asin(q*lamda/(4.0*pi)))

def theta2d(theta,lamda):
    return lamda/(2.0*sin(radians(theta/2.0)))

def d2theta(d,lamda):
    return 2.0*degrees(asin(lamda/(2.0*d)))

def q2d(q):
    return 2.0*pi/q

def d2q(d):
    return 2.0*pi/d


