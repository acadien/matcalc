#!/usr/bin/python
#Calculation tool for POSCARS and LAMMPS dumps.

import sys
#mine
import plotRemote as pr

from colors import float2rgb
from scipy import array
from numpy import linspace,savetxt
import matplotlib.cm as cm
import pylab as pl
from math import sqrt
import argparse
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


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description='''\
Apply some order parameter calculation to a POSCAR or LAMMPS Dump, then attempt to plot it or print out the result.\n
==========================================================
= Order Parameter must be one of:                        =
=   CN  : Coordination Number                            =
=   BO  : Bond Orientation (Q_l), requires -lval be set  =
=   RDF : Radial Distribution Function                   =
=   ADF : Angular Distribution Function                  =
=   BA  : Bond Angle Correlation, requires -gval be set  =
=   SF  : Structure Factor                               =
=   TET : Tetrahedral                                    =
=   TN  : Translational (tao)                            = 
==========================================================
''',epilog='''\
Examples:
charter.py RDF POSCAR1 POSCAR2 POSCAR2 -SA
charter.py ADF POSCAR -rcut 3.5
charter.py CN dumpfile.dat -N 25 -rcut 6.5
''')

#Need to implement test of file to determine if it is a lammps dump or poscar chart
#Also implement for OUTCAR files...

parser.add_argument('op',choices=orderParams.keys(),help="Order Parameter selection.")
parser.add_argument('fileNames',nargs='*',help="A list of POSCAR files (space delim)")
parser.add_argument('-rcut',help='Radius Cutoff distance when selecting neighbors, if not provided the first shell minimum is used.',type=float)
parser.add_argument('-lval',dest='lval',help='l-value, used in BO (Ql) calculation',type=int)
parser.add_argument('-gval',dest='lval',help='g-value, used in BA (g) calculation',type=int)
parser.add_argument('-S',dest='saveFlag',action='store_true',help='save to file instead of plotting',default='true')
parser.add_argument('-A',dest='averageFlag',action='store_true',help='average POSCAR results if possible')
parser.add_argument('-D',dest='debug',action='store_true',help='turn on debugging of neighbor selection')
#Need to implement this
parser.add_argument('-N',dest='lmpCfgReq',help='Configuration requested from LAMMPS Dump file',type=int)
#Need to implement this
parser.add_argument('-smooth',dest='windowSize',help='Turns on windowed averaging (smoothing) over all data sets',type=int,default=10)

args= parser.parse_args()
op = args.op
fileNames = args.fileNames
lval = args.lval

if op == "BA" and args.lval == None:
    print "Error: g-value (-gval) must be set for BA order parameter"
    exit(0)

if op == "BO" and args.lval == None:
    print "Error: l-value (-lval) must be set for BO order parameter"
    exit(0)

#Sort POSCAR names only based on the numbers in them
try:
    poscarNumbers=map(lambda x:float("".join(re.findall('\d+',x))),fileNames)
    if len(poscarNumbers) == len(fileNames):
        fileNames=zip(*sorted(zip(fileNames,poscarNumbers),key=lambda x:x[1]))[0]
except ValueError:
    pass

#Gather the order parameters
neighbs=None
orderVals=list()
for pn in fileNames:
    poscar=open(pn,"r").readlines()
    [basis,atypes,atoms,head,poscar] = poscarIO.read(poscar)

    #Get the order parameter
    orderVals.append(orderParams[op](array(atoms),array(basis),l=lval,neighbs=neighbs,rcut=args.rcut,debug=args.debug))

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

if args.averageFlag:
    if op not in ["BO","CN","TN","TET"]:
        t = zip(*[orderVals[i][1] for i in range(len(orderVals))])
        avgy = map(lambda x:sum(x)/len(x),t)
        avgx = orderVals[0][0]
    if args.saveFlag:
        if op in ["TN","TET"]:
            savetxt("AVERAGE."+op+str(lval),array([avgx,avgy]).T,delimiter=" ")
        elif op in ["ADF"]:
            savetxt("AVERAGE."+op+str(args.rcut),array([avgx,avgy]).T,delimiter=" ",header=" ".join(xylabels[op]),comments="")
        else:
            savetxt("AVERAGE."+op+str(lval),array([avgx,avgy]).T,delimiter=" ",header=" ".join(xylabels[op]),comments="")
        exit(0)
    else:
        pl.plot(avgx,avgy)

if args.saveFlag:
    for ov,pn in zip(orderVals,fileNames):
        if op in ["TN","TET"]:
            savetxt(pn+"."+op+str(lval),array(ov).T,delimiter=" ")

        elif op in ["ADF"]:
            if args.rcut==None:
                savetxt(pn+"."+op,array(ov).T,delimiter=" ",header=" ".join(xylabels[op]),comments="")
            else:
                savetxt(pn+"."+op+str(args.rcut),array(ov).T,delimiter=" ",header=" ".join(xylabels[op]),comments="")

        elif op in ["BO"]:
            vals,bins,dummy = pl.hist(ov,bins=int(sqrt(len(ov)))*5,normed=True,visible=False)
            savetxt(pn+"."+op+str(lval),array(zip(bins[:-1],vals)),delimiter=" ",header=" ".join(xylabels[op]),comments="")

        else:
            savetxt(pn+"."+op+str(lval),array(ov).T,delimiter=" ",header=" ".join(xylabels[op]),comments="")

    exit(0)

if op not in ["BO","CN","TN","TET"] and not args.averageFlag:
    for ov in orderVals:
        pl.plot(ov[0],ov[1])
    pl.xlabel(xylabels[op][0])
    pl.ylabel(xylabels[op][1])

if op=="BO":
    print "Average Bond Order <BO%d> ="%lval
    for i,ov in enumerate(orderVals):
        print fileNames[i],"\t\t",sum(ov)/len(ov)
        vals,bins,dummy = pl.hist(ov,bins=int(sqrt(len(ov)))*5,visible=False)#,histtype='step')
        pl.plot(array(bins[:-1]),array(vals),label=fileNames[i])
    pl.xlabel(xylabels["BO"][0])
    pl.ylabel(xylabels["BO"][1])

elif op=="CN":
    hd,rcuts=zip(*orderVals)
    pl.hist(hd,bins=range(0,16),normed=True,histtype='bar',align='left',rwidth=0.8)
    pl.xticks(range(min(map(min,hd)),max(map(max,hd))+1))
    for i,ov in enumerate(orderVals):
        cn=float(sum(ov[0]))/len(ov[0])
        print "Average CN (%s):"%fileNames[i],cn
        labels.append("CN%3.3f rcut%3.3f %s"%(cn,rcuts[i],fileNames[i]))
    pl.xlabel(xylabels["CN"][0])
    pl.ylabel(xylabels["CN"][1])

elif op=="TET":
    print "Average Tetrahedral Order <Sg> ="
    for i,ov in enumerate(orderVals):
        print fileNames[i],"\t\t",sum(ov)/len(ov)

elif op=="TN":
    print "Average Translational Order <tao> ="
    for i,ov in enumerate(orderVals):
        print fileNames[i],"\t\t",sum(ov)

if args.averageFlag:
    pl.legend(["Average "+" ".join(fileNames)],loc=0,fontsize=10)
else:
    if len(labels)>0:
        pl.legend(labels,loc=0,fontsize=10)
    else:
        pl.legend(fileNames,loc=0,fontsize=10)

#Don't plot for some values
if op not in ["TN","TET"]:
    pr.prshow("%s_chart_%s.png"%(sys.argv[1],op))
