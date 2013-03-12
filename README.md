This is code developed throughout my doctorate, developing semi-empirical models of materials and studying atomic structures in disordered materials. I'm sharing these programs with you to use as you please, feel free to modify, chop-up or ruin this code as much as you like.  Programs will tell you how to run them if run without any arguements, in general the name of the code should be self explanatory.

Required Libraries
==================
*   python-numpy
*   python-scipy
*   python-matplotlib

Optional Libraries
-------------------
*   python-gtk2
*   LAMMPS compiled with a python wrapper, with lammps.py in the PYTHONPATH

To Install
==========
Edit your bashrc (or equivalent) and add:
     matcalcdir="/path/to/matcalc/repo"
     matcalcsubdirs=${matcalcdir}:$(find ${matcalcdir} -type d | tr '\n' ':' | sed 's/:$//')
     export PYTHONPATH=$PYTHONPATH:${matcalcsubdirs}
     export PATH=$PATH:${matcalcsubdirs}

Description
===========
*   analysis - Generic tools for manipulating data and generating analyses with many applications e.g.: FFTs, band pass, multidimensional curve fitting, etc.
*   diffraction - Tools for x-ray pattern analysis and some curve fitting.
*   gmice - A couple of tools for running codes remotely on our HPC server.
*   moldyn - Scripts for reading lammps molecular dynamics output data, some analysis and plotting.
*   vasp - Scripts that plot and post-process VASP output data.
*   potfitting - Tools for potential fitting and analyzing potentials.
*   plotremote - tools for running these scripts from a remote server, makes working remotely easier when coupled with a file transfer system like dropbox.