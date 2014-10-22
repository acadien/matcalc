#!/usr/bin/python

import sys
#mine
import lammpsIO
import utils

utils.usage(["<input dump file>","<configuration #>","<output config file>"],3,3,"LAMMPS dump file must have bounds, atom type and atomic locations.")

iDump=sys.argv[1]
cfg=int(sys.argv[2])
bounds,types,atoms=lammpsIO.readConfig(iDump,cfg)
oCnfg=open(sys.argv[3],"w")
lammpsIO.dumpWriteConfig(oCnfg,bounds,types,atoms,"config made by %s, from file %s config %d"%(sys.argv[0].split("/")[-1],iDump,cfg))

