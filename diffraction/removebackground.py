#!/usr/bin/python
import sys
import pylab as pl
from scipy import interpolate
#mine
from chitools import *
from curvetools import *

#What does this script do?
#Given two chi files (input & background) takes the difference between the 
#two and dumps the result in the output chifile

#Arguement handline
def usage():
    print "%s <background chi file> <0:disable plot> <chifiles>"%sys.argv[0].split("/")[-1]

if len(sys.argv) < 4:
    usage()
    exit(0)

#Load arguements
bfilename=sys.argv[1]
enableplot=True
if sys.argv[2]=='0':
    enableplot=False

afilenames=sys.argv[3:]

for afilename in afilenames:
    ofilename="cleaned_"+afilename

    #Read in the achi file
    xxa,yya=readchi(afilename)
    xxb,yyb=readchi(bfilename)

    xxa2=[i for i in xxa]# if (i>xxb[0] and i<xxb[-1])]
    alow=xxa.index(xxa2[0])
    ahi=xxa.index(xxa2[-1])

    #Interpolate curve-b onto x-points of curve-a
    t=interpolate.splrep(xxb,yyb,k=5) #quintic spline
    yyb2=interpolate.splev(xxa2,t)

    #take the difference between the chifile and interpolated background
    ydiff=[yya[i+alow]-yyb2[i] for i in range(ahi-alow+1)]

    #Level out the difference by removing a linear portion
    #lcf=makefitlinear([xxa2[0],xxa2[-1]],[ydiff[0],ydiff[-1]])
    #leveldiff=[ydiff[i]-lineval(xxa2[i],lcf) for i in range(len(xxa2))]
    
    leveldiff=ydiff

    if enableplot:
        pl.figure()
        pl.plot(xxa,yya)
        pl.plot(xxa2,yyb2)
        pl.plot(xxa2,leveldiff)
        pl.plot([xxa[0],xxa[-1]],[0,0],ls='dotted')
        pl.title(afilename)
        pl.legend(["Input","Background","Difference"])
        pl.show()

    #Write the difference to the output file
    writechi(xxa2,leveldiff,ofilename)
    print "Successfully wrote to %s."%ofilename

