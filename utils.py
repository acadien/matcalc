#!/usr/bin/python

import sys

def usage(descrip,mn,mx,notes=None):
    nvars = len(sys.argv)-1
    if "-h" in sys.argv or nvars < mn or nvars > mx:
        print "%s %s"%(sys.argv[0].split("/")[-1]," ".join(descrip))
        if notes!=None:
            print notes
        exit(0)
