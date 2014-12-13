#!/usr/bin/python

import sys
#mine
from parserGens import parseLammpsColumns
import utils, lammpsIO

utils.usage(["lmpdump file"],1,1)

lmpFile = sys.argv[1]

nAtom = lammpsIO.nAtoms(lmpFile)
basis = lammpsIO.basis(lmpFile)

#find the header labels
for line in open(lmpFile):
    if "ITEM: ATOMS" in line:
        break
head = line.split()[2:]
print head

cfgIter = parseLammpsColumns(lmpFile,nAtom)

for i,config in enumerate(cfgIter):
    print i
