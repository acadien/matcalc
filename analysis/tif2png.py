#!/usr/bin/python

import sys
import pylab as pl
import Image
from math import sqrt

def usage():
    print "%s <tif image(s)>"%(sys.argv[0])

for fil in sys.argv[1:]:
    a=pl.array(Image.open(fil).getdata())
    
    a.shape=(2048,2048)

    avg=a.sum()/2048.0/2048.0

    imgplot=pl.imshow(a,vmin=0,vmax=10*avg)
    #pl.colorbar()
    pl.xticks([])
    pl.yticks([])
    imgplot.set_cmap('gist_earth')

    filout=fil.split(".")[0]+".png"    
    pl.savefig(filout)
