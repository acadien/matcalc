#!/usr/bin/python
import sys
import pylab as pl
import re
import os
import numpy as np
from numpy import *
#mine
from theta2q import theta2q

def smoothList(list,degree=5):  
    smoothed=np.zeros(len(list))  
    for i in range(len(smoothed)-(degree+1)/2):  
        smoothed[i+(degree+1)/2-1]=sum(list[i:i+degree])/float(degree)  
    return smoothed  

def localMax(seq,strt,delta):
    i = strt
    candidate = False
    m=-10000.0
    base=seq[strt]
    found=False
    for elem in seq[strt:]:  
        if elem-base >= delta: found=True
        if elem >= m:
            if not(candidate): base=elem
            m = elem        
            candidate = True
            i+=1
            continue
        else:      
            if found: return [i-1,m]
            m = elem        
            candidate = False 
            i+=1
            continue
    return [-1,m]    

def readchi(chifile):
    i=0
    for line in chifile:
        if i<3:
            i+=1
            continue
        line=line.lstrip()
        line=line.rstrip()
        if i==3:
            i+=1
            length=int(line)
            break
    xxs=np.zeros(length)
    yys=np.zeros(length)
    i=0;
    for line in chifile:
        line=line.strip().split()
        xxs[i]=line[0]
        yys[i]=line[1]
        i+=1
    return xxs,yys

#def theta2q(theta2,lamda):
    #return [4.0*pi*sin(radians(i/2.0))/lamda for i in theta2] #
#    return [lamda/2.0/sin(radians(i/2.0)) for i in theta2] 

##################
#Main Starts Here
##################
if len(sys.argv) < 2:
    print("Bad arguments, proper usage:")
    print("./plotchi.py <chifile1> <chifile2> ... <chifileN> <optional:lambda>")
    print "If the optional Lambda (Angstroms) arguement is given then the chifiles are plotted as a function of Q"
    exit(0)

try:
    lamda=float(sys.argv[-1])
    filenames=sys.argv[1:-1]
except ValueError:
    lamda=-1
    filenames=sys.argv[1:]

#####################
#Read in the chi file
xxs=list()
yys=list()
for filename in filenames:
    chifile=open(filename,'r')
    x,y=readchi(chifile)
    xxs.append(x)
    yys.append(y)

if lamda!=-1:
    print type(xxs[0])
    xxs=[map(lambda g:theta2q(g,lamda),xx) for xx in xxs]

"""
#####################
#Analyze the spectrum: find the peaks
if delta!=0:
    ind=0
    xxp=list()
    yyp=list()
#yysmth=yys
#yysmth=smoothList(yysmth)
    while True:
        [ind,pk]=localMax(yys,ind,delta)
        if ind==-1:
            break
        else:
            xxp.append(xxs[ind])
            yyp.append(yys[ind])
"""

#####################
#Plot
fig=pl.figure()
for x,y in zip(xxs,yys):
    pl.plot(x,y)
    #if delta!=0:
    #    pl.scatter(xxp,yyp)
if lamda!=-1:
    pl.xlabel("Q ($radians / \AA^{-1}$)")
else:
    pl.xlabel("2Theta (degrees)")
pl.ylabel("Intensity")
pl.legend(filenames)
pl.show()
        
    
