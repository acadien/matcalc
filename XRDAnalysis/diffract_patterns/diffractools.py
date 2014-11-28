#!/usr/bin/python
from math import *

#Returns the plane spacing given any lattice and lattice vectors
#Formulae taken from "X-Ray Diffraction" by Warren, page 21
#lattice constants: a,b,c
#lattice angles: A,B,C (given in degrees)
#lattice vector: string(hkl)
def calc_dspace(a,b,c,A,B,C,hkl):
    h,k,l=splithkl(hkl)
    A=radians(A)
    B=radians(B)
    C=radians(C)
    cA=cos(A)
    cB=cos(B)
    cC=cos(C)
    sA=sin(A)
    sB=sin(B)
    sC=sin(C)

    term1 = 1.0/(1.0 + 2.0*cA*cB*cC - cA*cA - cB*cB - cC*cC)
    term2 = h*h*sA*sA/a/a + k*k*sB*sB/b/b + l*l*sC*sC/c/c 
    term3 = 2.0*( h*k*(cA*cB-cC)/a/b + k*l*(cB*cC-cA)/b/c + l*h*(cC*cA-cB) )
    
    return sqrt(1.0/(term1*(term2+term3)))

#Given a wavelength and a d spacing, returns the expected theta value
def calc_theta(lamda,d):
    val=lamda/2.0/d
    if val>1.0 or val<-1.0:
        return -1
    return degrees(asin(val))

#Sorted by vector magnitude, from 1 to 35
hkls=["001","011","111","002","012","112","022","122","003","013","113","222","023","123","004","223","014","033","114","133","024","124","233","224","034","005","134","015","333","115","234","025","125","044","144","225","334","035","135"]

ortho_mults=[2,8,8,2,4,8,8,8,2,4,8,8,6,8,2,8,4,8,8,8,4,8,8,8,4,2,8,4,8,8,8,4,8,8,8,8,8,4,8]

def splithkl(hkl):
    return [int(i) for i in hkl]


"""          Multiplicities
Mag HKL      Cubic Ortho 
1   001      6     2        
2   011      12    8
3   111      8     8
4   002      6     2
5   012      24    4
6   112      24    8
7
8   022      12    8
9   122 003  30    10
10  013      24    4
11  113      24    8
12  222      8     8
13  023      24    4
14  123      48    8
15 
16  004      6     2
17  223 014  48    12
18  033 114  36    16
19  133      24    8
20  024      24    4
21  124      48    8
22  233      24    8
23 
24  224      24    8
25  034 005  30    6
26  134 015  72    12
27  333 115  32    16
28 
29  234 025  72    12
30  125      48    8
31 
32  044      12    8
33  144 225  48    16
34  334 035  48    12
35  135      48    8
"""

"""
System   Orthorhombic
hkl      8
hhl      8
hh0      8
0kk      8
hhh      8
hk0      4
h0l      4
0kl      4
h00      2
0k0      2
00l      2
Orthorhombic
8
8
8
8
8
4
4
4
2
2
2
Cubic
48
24
12
12
8
24
24
24
6
6
6
Tetragonal
16
8
4
8
8
8
8
8
4
4
2
Hexagonal
24
12
6
12
12
12
12
12
6
6
2
Monoclinic
4
4
4
4
4
4
2
4
2
2
2
Triclinic
2
2
2
2
2
2
2
2
2
2
2
"""
