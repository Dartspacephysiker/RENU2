#!/usr/bin/python

'''
Created on Jan 17, 2010

@author: wibble
'''
from struct import unpack
from sys import argv, exit, stderr, stdout

def help():
    stderr.write('''Usage: casc_extract.py <FI> <FS> <WI> <WS> <IN> [OUT]
    FI : Frame Interval, 1 if word appears in every minor frame.
    FS : Frame Start, first minor frame word appears in.
    WI : Word Interval, 75 if word appears once per minor frame.
    WS : Word Start, where in the minor frame this word first appears.
    IN : Input file
    OUT : Output file [optional, defaults to stdout].

''')

if len(argv) < 6 or '-h' in argv or '--help' in argv or '/h' in argv:
    help()
    exit()

wpf = 75 # words per frame
bpff = 89 # bytes per file frame
tpw = 0.0000016666666667 # time per word

frame_int = int(argv[1])
frame_start = int(argv[2])
word_int = int(argv[3])
word_start = int(argv[4])
stderr.write('Input: Binary from ' + argv[5] + '.\n')
infile = open(argv[5], 'rb')

if len(argv) < 7:
    ofile = stdout
    stderr.write('Output: ASCII to standard output.\n')
else:
    ofile = open(argv[6] + ".asc", 'w+')
    tfile = open(argv[6] + ".times", 'w+')
    gfile = open(argv[6] + ".gps", 'w+')
    stderr.write('Output: ASCII to ' + argv[6] + '.\n')

def bstr(n):
	return ''.join([str(tframe[15] >> x & 1) for x in (7,6,5,4,3,2,1,0)])

minfr = 0
minfrl = 0
gpsbyte = 0
gpsbits = 0
while True:
    frame = infile.read(bpff)
    
    if len(frame) != bpff:
        break
    
    header = unpack('>HHHH', frame[0:8])
    
    if header == (12291, 816, 1, 20736):  # 3003 0330 0001 5100
        # File synced
        minfr += 1
        
        tstamp = unpack('>Q', "\x00\x00" + frame[8:14])[0]
        tframe = [x for x in unpack('>' + str(wpf) + 'B', frame[14:])]
        
        sec = ((tstamp >> 10) & 134217727) / 1000.0 # 27-bit milliseconds-of-day
        sec += (tstamp & 1023) * 0.000001 # 10-bit microseconds-of-millisecond
        sec -= 39840 # Subtract launch time to get seconds-since-launch (11:04:00 UT)
        sec += (word_start - 1) * tpw
        
        if not tframe[0:3] == [0xFA, 0xF3, 0x20]:
            print "Missing frame sync at frame ", minfr, "."
            continue
        
        if not (tframe[15] & 0b11100000 == 128 or tframe[15] & 0b11100000 == 192):
            print "Bad SFID at frame {0}: {1}, byte {2}.".format(minfr,bstr(tframe[15]),infile.tell())
            continue
        
        nframe = (tframe[15] & 0b00011111) + 1 # extract frame counter
        if (nframe - frame_start) % frame_int == 0:
            # correct minor frame
            for j in range(word_start - 1 + 3, wpf, word_int):
                mword = tframe[j]
                lword = tframe[j + 1]
                
                gpsbyte = (gpsbyte<<1)|(lword&0b00000001)
                gpsbits += 1
                if gpsbits == 8:
                    gfile.write("{0}\n".format(bstr(gpsbyte)))
                    gpsbyte = 0
                    gpsbits = 0
#                if lword&0b00000100 != 4:
#                    print "Bad data (bit 5) at word {0}, frame {1}.".format(j, minfr)
#                    continue
                
                word = (mword<<4)|(lword&0b11110000)
                
                tfile.write("{0:.7f}\n".format(sec))
                ofile.write("{0:-4d}\n".format(word))
                if word_int != wpf:
                    sec += word_int * tpw
        
    else:
        stderr.write("File sync lost at byte " + str(infile.tell()) + ".\n")
        exit()
    
   # frame = infile.read(bpff)

ofile.close()
