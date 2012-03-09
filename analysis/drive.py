#!/usr/bin/python

import itertools,random
from math import *
from numpy import array
#mine
from struct_tools import *

random.seed(3)

atoms=itertools.product(range(8),repeat=3)
atoms=[[(j+random.random()*3.0-1.5)%10.0 for j in i] for i in atoms]
atoms=[[(j+random.random()*3.0-1.5)%10.0 for j in i] for i in atoms]
atoms=[[(j+random.random()*3.0-1.5)%10.0 for j in i] for i in atoms]

atoms.append([0,0,0])
atoms.append([0,0,1])
atoms.append([0,1,0])
atoms.append([0,1,1])
atoms.append([1,0,0])
atoms.append([1,0,1])
atoms.append([1,1,0])
atoms.append([1,1,1])
atoms=array([array(i) for i in atoms])
N=len(atoms)

neighbs = neighbors(atoms,2.0)
squares = find_squares(atoms,neighbs,atol=1) 
cubes = find_cubes(squares)
print squares
print cubes

