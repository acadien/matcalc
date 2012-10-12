#!/bin/bash

#This script adds all of the configurations in the finals directory to the force database.

dirs=`ls ./finals/`

for d in $dirs
do
    ./append_ocdb.py $1 ./finals/$d 1 1
done
