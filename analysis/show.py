#!/usr/bin/python

import sys,Image
import pylab as pl
from numpy import *

def usage():
    print "%s <tif/png/jpeg/eps image> <optional:min,max>"%(sys.argv[0])


def loadimg(name):
    picdata=Image.open(name)
    a,b,c,d=picdata.getbbox()
    width=c-a
    height=d-b+1
    return array(picdata.getdata()).reshape([width,height])

if len(sys.argv)==0:
    usage()
    exit(0)

picdata=loadimg(sys.argv[1])
if len(sys.argv)==3:
    a,b=map(float,sys.argv[2].split(","))
    imgplot=pl.imshow(picdata,vmin=a,vmax=b)
    imgplot.set_cmap('gist_earth')
else:
    imgplot=pl.imshow(picdata)

pl.show()

