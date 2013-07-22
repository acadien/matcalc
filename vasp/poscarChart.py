#!/usr/bin/python


import sys
#mine
import plotRemote as pr

from colors import float2rgb
from scipy import array
from numpy import linspace,savetxt
import matplotlib.cm as cm
import pylab as pl
from math import sqrt
import re
#mine

from orderParam import *
import poscarIO

orderParams={"CN":  coordinationNumber, 
             "BO":  bondOrientation , 
             "RDF": radialDistribution , 
             "ADF": angleDistribution , 
             "BA":  bondAngleCorr , 
             "SF":  structureFactor, 
             "TET": tetrahedral, 
             "TN":  translational}

def usage():
    print "%s <order parameter> <POSCAR Files (space delim)> <S = save to file instead of plotting> <A = average POSCAR results if possible>"%sys.argv[0].split("/")[-1]
    print "Order Parameter must be one of:"
    print "   CN#.#  : Coordination Number with rcutoff (optional)"
    print "   BO# : Bond Orientation (Q) with l=#"
    print "   RDF : Radial Distribution Function"
    print "   ADF : Angular Distribution Function"
    print "   BA# : Bond Angle Correlation (g)"
    print "   SF  : Structure Factor"
    print "   TET : Tetrahedral"
    print "   TN  : Translational (tao)"
    print ""

if len(sys.argv) < 3:
    usage()
    exit(0)

lval=None
op = sys.argv[1]
if op[:2] in ["BO","BA"]:
    lval=int(op[-1])
    op=op[:2]

rcut=None
if op[:2] == "CN" and len(op)>2:
    rcut=float(op[2:])
    op=op[:2]

orderVals=list()
poscarNames=sys.argv[2:]
saveFlag=False
averageFlag=False

if poscarNames[-1] in ["SA","AS"]:
    poscarNames.pop()
    saveFlag=True
    averageFlag=True
elif poscarNames[-1]=="S":
    poscarNames.pop()
    saveFlag=True
elif poscarNames[-1]=="A":
    poscarNames.pop()
    averageFlag=True

#Sort POSCAR names only based on the numbers in them
try:
    poscarNumbers=map(lambda x:float("".join(re.findall('\d+',x))),poscarNames)
    if len(poscarNumbers) == len(poscarNames):
        poscarNames=zip(*sorted(zip(poscarNames,poscarNumbers),key=lambda x:x[1]))[0]
except ValueError:
    pass

#Gather the order parameters
l=lval
neighbs=None
debug=False
for pn in poscarNames:
    poscar=open(pn,"r").readlines()
    [basis,atypes,atoms,head,poscar] = poscarIO.read(poscar)

    #Get the order parameter
    orderVals.append(orderParams[op](array(atoms),array(basis),l=l,neighbs=neighbs,rcut=rcut,debug=debug))

#======================================================
#                       Plot!
#======================================================
labels=[]
if lval==None:
    lval=""
xylabels={"BO" :[r"BondOrder($Q_%s$)"%str(lval), r"$P(Q_%s)$"%str(lval)],
          "CN" :[r"Coordination_Number",         r"P(CN)"],
          "RDF":[r"R$(\AA)$",                    "#_Bonds"],
          "ADF":[r"$\theta(deg)$",               "#_Bond_Angles"],
          "BA" :[r"R$(\AA)$",                    r"$G_{%s}(r)$"%str(lval)],
          "SF" :[r"Q$(\AA^{-1})$",               r"S(Q)"]}

if averageFlag:
    if op not in ["BO","CN","TN","TET"]:
        t = zip(*[orderVals[i][1] for i in range(len(orderVals))])
        avgy = map(lambda x:sum(x)/len(x),t)
        avgx = orderVals[0][0]
    if saveFlag:
        if op in ["TN","TET"]:
            savetxt("AVERAGE."+op+str(lval),array([avgx,avgy]).T,delimiter=" ")
        else:
            savetxt("AVERAGE."+op+str(lval),array([avgx,avgy]).T,delimiter=" ",header=" ".join(xylabels[op]),comments="")
        exit(0)
    else:
        pl.plot(avgx,avgy)

if saveFlag:
    for ov,pn in zip(orderVals,poscarNames):
        if op in ["TN","TET"]:
            savetxt(pn+"."+op+str(lval),array(ov).T,delimiter=" ")
        else:
            savetxt(pn+"."+op+str(lval),array(ov).T,delimiter=" ",header=" ".join(xylabels[op]),comments="")
    exit(0)

if op not in ["BO","CN","TN","TET"] and not averageFlag:
    for ov in orderVals:
        pl.plot(ov[0],ov[1])
    pl.xlabel(xylabels[op][0])
    pl.ylabel(xylabels[op][1])

if op=="BO":
    print "Average Bond Order <BO%d> ="%lval
    for i,ov in enumerate(orderVals):
        print poscarNames[i],"\t\t",sum(ov)/len(ov)
        vals,bins,dummy = pl.hist(ov,bins=int(sqrt(len(ov)))*5,normed=True,visible=False)#,histtype='step')
        pl.plot(array(bins[:-1]),array(vals),label=poscarNames[i])
    pl.xlabel(xylabels["BO"][0])
    pl.ylabel(xylabels["BO"][1])

elif op=="CN":
    hd,rcuts=zip(*orderVals)
    pl.hist(hd,bins=range(0,16),normed=True,histtype='bar',align='left',rwidth=0.8)
    pl.xticks(range(min(map(min,hd)),max(map(max,hd))+1))
    for i,ov in enumerate(orderVals):
        cn=float(sum(ov[0]))/len(ov[0])
        print "Average CN (%s):"%poscarNames[i],cn
        labels.append("CN%3.3f rcut%3.3f %s"%(cn,rcuts[i],poscarNames[i]))
    pl.xlabel(xylabels["CN"][0])
    pl.ylabel(xylabels["CN"][1])

elif op=="TET":
    print "Average Tetrahedral Order <Sg> ="
    for i,ov in enumerate(orderVals):
        print poscarNames[i],"\t\t",sum(ov)/len(ov)

elif op=="TN":
    print "Average Translational Order <tao> ="
    for i,ov in enumerate(orderVals):
        print poscarNames[i],"\t\t",sum(ov)

if averageFlag:
    pl.legend(["Average "+" ".join(poscarNames)])
else:
    if len(labels)>0:
        pl.legend(labels,loc=0)
    else:
        pl.legend(poscarNames,loc=0)
#Don't plot for some values
if op not in ["TN","TET"]:
    pr.prshow("%s_chart_%s.png"%(sys.argv[1],op))
