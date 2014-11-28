#!/usr/bin/python

import os
import re
from ID3 import *

path='./'
dirlist=os.listdir(path)
pardir=os.getcwd()
album=re.search(r'(/[\w,\s]+)+',pardir)
album=album.group(1).lstrip('/')
for fname in dirlist:
    try:
        m = re.search('(\d\d)-(\w+)-(\w+)_www\.file',fname)
        id3info = ID3(fname)
    except:
        continue
    print id3info
    id3info['TRACKNUMBER'] = m.group(1)
    artist = m.group(2)
    id3info['ARTIST'] = re.sub('_',' ',artist).capitalize()
    song = m.group(3)
    id3info['SONG']=re.sub('_',' ',song).capitalize()
    id3info['ALBUM']=album
    #print track+artist+song
    #convert='mp3info -f -t '+song+' -n '+track+' -a '+artist+' -l '+album+' '+fname
    #os.system(convert)
