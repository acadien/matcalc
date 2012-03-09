#!/usr/bin/python

from math import *

def theta2q(2theta,lamda):
    return [4.0*pi*sin(i/2.0)/lamda for i in 2theta] 
