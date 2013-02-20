#!/usr/bin/python

import sys

#bandE,dos,tdos = read(doscar)
#if spin is enabled then dos and tdos will have spin up and spin down components
def read(doscar):
    Natoms=int(doscar[0].split()[0])
    efermi = float(doscar[5].split()[3])
    doscar=doscar[5:]
    NDOS=int(doscar.pop(0).split()[2])
    data=zip(*[map(float,i.split()) for i in doscar[1:NDOS]])
    bandE=data[0]
    if len(data)==3:
        dos=data[1]
        tdos=data[2]
    elif len(data)==5:
        dos=data[1:3]
        tdos=data[3:5]
    return efermi,bandE,dos,tdos
