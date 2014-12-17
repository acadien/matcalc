#!/usr/bin/python

import sys
import operator
#mine
from parserGens import parseLammpsColumns
import utils, lammpsIO

#Given a lammps dump and some criteria, dumps the atoms that have valid criteria
#also gives easy access to their coordinates for further analysis

utils.usage(["<lmpdump file>","<selection crit>","lt/gt val eg: gt0.8"],3,3)

lmpFile = sys.argv[1]
selHead = sys.argv[2]
criteria = sys.argv[3]

outFileAtoms = lmpFile+"_"+selHead+"_"+criteria+".atoms"
outFileSize = lmpFile+"_"+selHead+"_"+criteria+".len"

op = operator.lt
if "gt" in criteria:
    op = operator.gt
val = float(criteria[2:])
def condition(compareVal):
    return op(compareVal,val)

nAtom = lammpsIO.nAtoms(lmpFile)
basis = lammpsIO.basis(lmpFile)

#find the header labels
for line in open(lmpFile):
    if "ITEM: ATOMS" in line:
        break
head = line.split()[2:]

if selHead not in head:
    print "Selection Criteria is not available, choose from:"
    print head
    exit(0)

selIndex = head.index(selHead)
cfgIter = parseLammpsColumns(lmpFile,nAtom)

print "Writing %s and %s"%(outFileAtoms,outFileSize)
outAtoms = open(outFileAtoms,"w")
outSize = open(outFileSize,"w")
for i,config in enumerate(cfgIter):
    try:
        criteriaColumn = zip(*config)[selIndex]
    except IndexError:
        break

    #The full cluster with all of its extraneous information... oh boy
    clusterFull = [a for j,a in enumerate(config) if condition(criteriaColumn[j])]


    cluster = [int(a[0]) for a in clusterFull]
    outAtoms.write(" ".join(map(str,cluster))+"\n")
    outSize.write(str(len(cluster))+"\n")
    
