#!/usr/bin/python

#Takes the contents from the N'th "ITEM: TIMESTEP" from a lammps dump and sticks it in a new file

import sys

def usage():
    print "%s <lammps dump with atom cfg> <desired cfg (-1 for last)>"%sys.argv[0].split("/")[-1]

if len(sys.argv)<2:
    usage()
    exit(0)
    
cfg=int(sys.argv[2])

cnt=0
foundFlag=False
fp=open(sys.argv[1],"rb")
copydata=list()

if cfg==-1:
    for line in fp:
        if line=="ITEM: TIMESTEP\n":
            cfg+=1
    fp.seek(0)

for i,line in enumerate(fp):
    if line=="ITEM: TIMESTEP\n":
        if cnt==cfg:
            foundFlag=True
            copydata.append(line)
            break
        cnt+=1
        continue

if not foundFlag:
    print "PARSE ERROR"
    print "Requested configuration not found, either too large or not valid int."
    exit(0)


for line in fp:
    if line=="ITEM: TIMESTEP\n":
        break
    copydata.append(line)
open("cfg%4.4d_"%cnt+sys.argv[1],"w").writelines(copydata)




