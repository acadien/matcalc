#!/usr/bin/python

import sys
import pylab as pl
from matplotlib import image

def usage():
    print "%s <tif/png/jpeg/eps image> <optional: 0-grayscale, 1-color>"%(sys.argv[0])

picdata=image.imread(sys.argv[1])
if len(sys.argv)==3 and sys.argv[2]=="1": 
    imgplot=pl.imshow(picdata[:,:,0])
    imgplot.set_cmap('spectral')
else:
    imgplot=pl.imshow(picdata)

pl.show()
