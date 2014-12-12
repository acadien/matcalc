#!/usr/bin/python

import sys
#mine
import utils
import parserGens
import datatools

utils.usage(["<.tetraFile>","<window size (int or all)>"],2,2)

inputFile = sys.argv[1]
if ".tetra" not in inputFile:
    print "wrong input file dummy."
    exit(0)

windowSize = sys.argv[2]
outputFile = inputFile+".avg"+windowSize

configIterator = parserGens.parseEnsemble(inputFile,keepAvg=True)
tetras = zip(*[a for a in configIterator])

if windowSize == "all":
    tetraAvg = [sum(atet)/len(atet) for atet in tetras]
    tetraLines = [" ".join(map(str,tetraAvg))+"\n"]
else:
    windowSize = int(windowSize)
    tetraAvg = zip(*[datatools.windowAvg(atet,n=windowSize,option='valid') for atet in tetras])
    tetraLines = [" ".join(map(str,line))+"\n" for line in tetraAvg]

out = open(outputFile,"w")
out.write("AvgTetra TetraPerAtom\n")
out.writelines(tetraLines)
