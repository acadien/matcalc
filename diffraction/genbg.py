#!/usr/bin/python
import sys
import pylab as pl
from numpy import *
from scipy import interpolate
from math import *
#mine
from curvetools import *
from chitools import *

def stddev(items):
    m=sum(items)/len(items)
    return sqrt(sum([fabs(i-m)**2 for i in items])/len(items))

#Finds a bunch of minima in the range defined.
#Works using a paved minimum method
def findmins(yys,start,end,number):
    temp=[i for i in yys[:]]
    m=int((end-start)/number)-5
    points=list()
    points_i=list()
    while len(points)<number:
        p=min(temp[start:end])
        p_i=temp.index(p)
        points.append(p)
        points_i.append(p_i)
        mn=max(start,p_i-m)
        mx=min(end,p_i+m+1)
        temp[mn:mx]=[temp[i]+50 for i in range(mn,mx)]
    return points_i

#################
#Helper functions
#################
funcval={"line":lineval,
         "quad":quadval,
         "gaus":gaussianval}

#Finds outliers relative to some fitted curve given
def findoutliers(xs,ys,ctyp,args):
    #Find the standard deviation from the function given
    n=len(xs)
    sd=sum([fabs(ys[i]-funcval[ctyp](xs[i],args))**2 for i in range(n)])/n
    return [i for i in range(n) if fabs(ys[i]-funcval[ctyp](xs[i],args))>sd]

#does a least squares fit on the desired function
def dofit(bgx,bgy,ctyp):
    if ctyp=="line":
        return makefitlinear(bgx,bgy)
    elif ctyp=="quad":
        return makefitquadratic(bgx,bgy)
    elif ctyp=="gaus":
        half=len(bgx)/2           
        return makefitgaussian(bgx,bgy,(1,bgx[half],1))

def genbg(xs,cf,ctyp):
    if ctyp=="line":
        return [lineval(x,cf) for x in xs]
    elif ctyp=="quad":
        return [quadval(x,cf) for x in xs]
    elif ctyp=="gaus":
        return [gaussianval(x,cf) for x in xs]

def usage():
    print "%s <ichifile> <ochifile> <start2theta> <end2theta> <line|quad|gaus> <specific points to keep (from previous run)>"%sys.argv[0]
    print "Given input chifile, start/stop locations, generates the background curve requested (linear/quadratic/gaussian) and writes it to the output chi file"

#####
#Main
if len(sys.argv) < 6:
    usage()
    exit()

ifilename=sys.argv[1]
ofilename=sys.argv[2]
s2t=float(sys.argv[3])
e2t=float(sys.argv[4])
ctyp=sys.argv[5]
tokeep=list()
if len(sys.argv)>6: 
    tokeep=[int(i) for i in sys.argv[6:]]

#Read in the chi file
xxs,yys=readchi(ifilename)
    
#Get the starting and ending 2theta indeces
s2ti= sum([1 for i in xxs if i<s2t]) 
if e2t<0:
    e2ti= len(xxs)-1
else:
    e2ti= sum([1 for i in xxs if i<e2t])

#Remove a linear fit from the background.
lcf=makefitlinear([xxs[s2ti],xxs[e2ti]],[yys[s2ti],yys[e2ti]])
levelxxs=xxs[s2ti:e2ti]
levelyys=[yys[i]-lineval(xxs[i],lcf) for i in range(s2ti,e2ti)]

#Find a bunch of minima scattered throughout the curve
startpoints=10
bgi=list(set(findmins(levelyys,0,len(levelyys)-1,startpoints)))
bgi=sorted(bgi)
bgx=[levelxxs[i] for i in bgi]
bgy=[levelyys[i] for i in bgi]

#Do the initial fit and eliminate outliers, gradually improve the fit
"""
for i in range(5):
    coef=dofit(bgx,bgy,ctyp)
    outliers=findoutliers(bgx,bgy,ctyp,coef)
    outliers=reversed(sorted(outliers)) #have to sort/reverse to ensure pop works
    #Remove outliers that aren't local minima
    for out in outliers:
        i=bgi[out]
        y=levelyys[i]
        mn=max(0,i-25)
        mx=min(len(levelyys)-1,i+25)
        if not(y==min(levelyys[mn:mx])):
            bgi.pop(out) 
    bgx=[levelxxs[i] for i in bgi]
    bgy=[levelyys[i] for i in bgi]
    if len(bgi)<10:
        break
"""

if len(tokeep)>0:
    bgi=[bgi[i] for i in range(len(bgi)) if i in tokeep]
#    bgi.append(max([j for j,i in enumerate(levelxxs) if i > 15.33 and i < 15.45]))
    bgx=[levelxxs[i] for i in bgi]
    bgy=[levelyys[i] for i in bgi]

#Do the final fit
cf=dofit(bgx,bgy,ctyp)

#levelxxs=xxs[s2ti:e2ti-15]
shift=0
totalbackground=[funcval[ctyp](x,cf)+lineval(x,lcf) for x in levelxxs]
levelbackground=[funcval[ctyp](x,cf)-shift for x in levelxxs]

#writechi(levelxxs,levelbackground,ofilename)
writechi(levelxxs,totalbackground,ofilename)
print "The curve shown here was written to %s"%ofilename
    
#Plot for testing
pl.figure()
pl.plot(xxs,yys)
[pl.scatter(levelxxs[i],totalbackground[i],c="yellow") for i in bgi]
pl.plot(levelxxs,totalbackground)
[pl.text(levelxxs[i],totalbackground[i]-20,str(j),fontsize=10) for j,i in enumerate(bgi)]
pl.show()



