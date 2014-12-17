#!/usr/bin/python

import sys
#mine
import utils
import parserGens
import datatools

utils.usage(["<.ensembleFile>","<window size (int or all)>"],2,2)

inputFile = sys.argv[1]
possibleSuffixes = [".tetra",".rmsd",".cn"]
if not any([suffix for suffix in possibleSuffixes if suffix in inputFile]):
    print "wrong input file dummy."
    exit(0)

windowSize = sys.argv[2]
outputFile = inputFile+".avg"+windowSize

#assume first line is a header
head = open(inputFile,"r").readline()

configIterator = parserGens.parseEnsemble(inputFile,keepAvg=True)
ops = zip(*[a for a in configIterator])

if windowSize == "all":
    opAvg = [sum(op)/len(op) for op in ops]
    opLines = [" ".join(map(str,opAvg))+"\n"]
else:
    windowSize = int(windowSize)

    opAvg = zip(*[datatools.wsmooth(atet,windowSize) for atet in ops])
    opLines = [" ".join(map(str,line))+"\n" for line in opAvg]

out = open(outputFile,"w")
out.write(head)
out.writelines(opLines)
