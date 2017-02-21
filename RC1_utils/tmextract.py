#!/usr/bin/python

'''
Created on Jan 17, 2010

@author: wibble
'''
from struct import unpack
from sys import argv,exit,stderr,stdout

def help():
    stderr.write('''Usage: tmextract.py <FI> <FS> <WI> <WS> <IN> [OUT]
    FI : Frame Interval, 1 if word appears in every minor frame.
    FS : Frame Start, first minor frame word appears in.
    WI : Word Interval, 160 if word appears once per minor frame.
    WS : Word Start, where in the minor frame this word first appears.
    IN : Input file
    OUT : Output file [optional, defaults to stdout].

''')

if len(argv) < 6 or '-h' in argv or '--help' in argv or '/h' in argv:
    help()
    exit()

frame_int = int(argv[1])
frame_start = int(argv[2])
word_int = int(argv[3])
word_start = int(argv[4])
stderr.write('Input: Binary from ' + argv[5] + '.\n')
infile = open(argv[5],'rb')

if len(argv) < 7:
    outfile = stdout
    stderr.write('Output: ASCII to standard output.\n')
else:
    outfile = open(argv[6],'w+')
    stderr.write('Output: ASCII to ' + argv[6] + '.\n')

frame = infile.read(334)

minfr = 0
minfrl = 0
while len(frame) == 334:
    header = unpack('>HHHH',frame[0:8])
    
    if header == (12291,816,1,17921):  # 3003 0330 0001 4601
        # File synced
        minfr += 1

        tstamp = unpack('>Q',"\x00\x00"+frame[8:14])[0]
        tframe = [x&1023 for x in unpack('>'+'H'*160,frame[14:])]
        
        sec = ((tstamp>>10)&134217727)/ 1000.0 # 27-bit milliseconds-of-day
        sec += (tstamp&1023)*0.000001 # 10-bit microseconds-of-millisecond
        sec -= 35350.234 # Subtract launch time to get seconds-since-launch

        if tframe[0:2] == [1003,819]:
            # synced to Frame
#            print "frsync",frame[30:47] 
            if minfrl != 0:
                stderr.write("Lost minor frame sync at {0} ({1} lost).\n".format(minfrl,minfr-minfrl))
                minfrl = 0

        nframe = (tframe[3] & 31) + 1 # extract frame counter
        if (nframe-frame_start) % frame_int == 0:
            for j in range(word_start-1+3,160,word_int):
                word = tframe[j]
                # we have 10-bit words, and for some reason the PTP format does not
                # fill unused bits with zeroes, so we need to mask those bits out 
                
                outfile.write("{0:f} {1:-4d}\n".format(sec,word))
                if word_int != 160:
                    sec += word_int*0.00000125
                
    else:
        stderr.write("File sync lost at byte " + str(infile.tell()) + ".\n")
        exit()

    frame = infile.read(334)
    
outfile.close()