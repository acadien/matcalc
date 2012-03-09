#!/usr/bin/python

import os

def dir(localdir,remotedir):
    os.system("scp -r %s gmice.gmu.edu:/home/acadien/%s/"%(localdir,remotedir))
