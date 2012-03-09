#!/usr/bin/python

#This script takes a completed vasp run and generates a new directory with the same run generating a POSCAR from the last stage of the previous OUTCAR.

import sys,os
#mine
import outcar2poscar

def usage():
    print "%s <old_dir> <new_dir>"%sys.argv[0]
    exit(0)

if len(sys.argv)!=3:
    usage()

stage1=sys.argv[1].strip("/")+"/"
stage2=sys.argv[2].strip("/")+"/"

os.system("mkdir %s"%stage2)
os.system("cp %sKPOINTS %s"%(stage1,stage2))
os.system("cp %sINCAR %s"%(stage1,stage2))
os.system("cp %sPOTCAR %s"%(stage1,stage2))
#os.system("cd %s"%stage1)

outcar2poscar.outcar2poscar(stage1,"%sPOSCAR"%stage2,-1)

