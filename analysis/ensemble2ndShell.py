#!/usr/bin/python

import sys
import itertools
#mine
import utils,parserGens,datatools

utils.usage(["<file.neighbs>","<ensemble file (tetra, cn, rmsd, w/e)>"],2,2)

neighbFile = sys.argv[1]
ensembFile = sys.argv[2]

neighbGen = parserGens.parseNeighbor(neighbFile)
ensembGen = parserGens.parseEnsemble(ensembFile)

#secondShellEnsemb=list()
efp = ensembFile.split(".")
outfile = ".".join(efp[0:-1]) + ".2shell" + efp[-1]
out = open(outfile,"w")
out.write("2ndShellAvg%s 2ndShellPerAtom%s\n"%(efp[1],efp[1]))
stopFlag=False
for ensembG,neighbG in itertools.izip(ensembGen,neighbGen):

    secondShellEnsemb = list()
    for a,firstNeighbs in enumerate(neighbG):

        secondNeighbs = [neighbG[i] for i in firstNeighbs]
        totalNeighb = set(datatools.flatten([firstNeighbs]+secondNeighbs))
        N = len(totalNeighb)        
        if N==0:
            secondShellEnsemb.append(ensembG[a])
        else:
            try:
                secondShellEnsemb.append(sum([ensembG[i] for i in totalNeighb])/N)
            except IndexError:
                stopFlag=True
                break
    if stopFlag: break
    avg = sum(secondShellEnsemb)/len(secondShellEnsemb)    
    out.write(str(avg)+" " + " ".join(map(str,secondShellEnsemb))+"\n")

out.close()
print "Wrote %s."%outfile

