#!/usr/bin/python
import sys,os

class _Getch:
    """Gets a single character from standard input.  Does not echo to the
screen."""
    def __init__(self):
        self.impl = _GetchUnix()

    def __call__(self): return self.impl()

class _GetchUnix:
    def __init__(self):
        import tty, sys

    def __call__(self):
        import sys, tty, termios
        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        try:
            tty.setraw(sys.stdin.fileno())
            ch = sys.stdin.read(1)
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
        return ch

dirs=["./"]
if len(sys.argv)>1:
    dirs=[i.strip("/")+"/" for i in sys.argv[1:]]

getch = _GetchUnix()

print "Are you sure you want to delete all VASP output files (keeping inputs)?[N/Y]:"
print dirs
response= getch()

if response == "Y":
    files=["AECCAR0","AECCAR1","AECCAR2","CHG","CHGCAR","CONTCAR","DOSCAR","EIGENVAL","ELFCAR","IBZKPT","OSZICAR","OUTCAR","output","PCDAT","PROCAR","vasprun.xml","WAVECAR","XDATCAR"]
    for d in dirs:
        for f in files:
            try:
                os.remove(d+f)
            except OSError:
                pass
else:
    print "Exiting."
