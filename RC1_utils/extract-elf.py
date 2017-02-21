#!/usr/bin/env python

'''
Created on Jan 17, 2010

@author: wibble
'''
from struct import unpack
from sys import argv,exit,stderr,stdout

def help():
	stderr.write('''Usage: tmextract.py <IN> [OUT]
	IN : Input file
	OUT : Output file [optional, defaults to stdout].

''')

if len(argv) < 2 or '-h' in argv or '--help' in argv or '/h' in argv:
	help()
	exit()

stderr.write('Input: Binary from ' + argv[1] + '.\n')
infile = open(argv[1],'rb')

outfile = open(argv[2],'w+')
stderr.write('Output: ASCII to ' + argv[2] + '.\n')
	
outfile.write("# Time\t\tELF_x_lo\tELF_x_hi\tELF_y_lo\tELF_y_hi\tELF_s_lo\tELF_s_hi\n")

frame = infile.read(100000)

elfstruct = [2**17]*6 # samples are 10-bit, so 2^17 is OOB/error
minfr = 0
minfrl = 0
filesync = "\x30\x03\x03\x30\x00\x01\x46\x01"

cpos = frame.find(filesync)
frame = frame[cpos:]

while len(frame) >= 334:
	# File synced
	minfr += 1

	tstamp = unpack('>Q',"\x00\x00"+frame[8:14])[0]
	tframe = [x&1023 for x in unpack('>'+'H'*160,frame[14:334])]
	
	sec = ((tstamp>>10)&134217727)/ 1000.0 # 27-bit milliseconds-of-day
	sec += (tstamp&1023)*0.000001 # 10-bit microseconds-of-millisecond
	sec -= 35350.234 # Subtract launch time to get seconds-since-launch

	if tframe[0:2] == [1003,819]:
		# synced to Frame
#			 print "frsync",frame[30:47] 
		if minfrl != 0:
			stderr.write("Lost minor frame sync at {0} ({1} lost).\n".format(minfrl,minfr-minfrl))
			minfrl = 0

	nframe = (tframe[3] & 31)  # extract frame counter
	subfr = nframe%5
	
#	print 'a',subfr, elfstruct
	if subfr == 2: # new timestamp -- check, output, and reset
		if max(elfstruct) != 2**17: # struct looks filled
#			print sec, elfstruct
			outfile.write(("{:f}"+"\t{:-4d}\t"*6+"\n").format(prevsec, *elfstruct)) # print cluster
		else:
			print sec, elfstruct # output error times and values
			
		elfstruct = [2**17]*6 # reset
#		print elfstruct
		prevsec = sec
		
		# [xlo, xhi, ylo, yhi, slo, shi]
	if subfr == 0:
		elfstruct[2] = tframe[119-1+3]
	elif subfr == 1:
		elfstruct[3] = tframe[119-1+3]
	elif subfr == 2:
		elfstruct[4] = tframe[119-1+3]
	elif subfr == 3:
		elfstruct[0] = tframe[118-1+3]
		elfstruct[5] = tframe[119-1+3]
	elif subfr == 4:
		elfstruct[1] = tframe[118-1+3]

#	print subfr, elfstruct
#				 if word_int != 160:
#					 sec += word_int*0.00000125

	cpos = frame.find(filesync, 1)
#	print cpos
	
	if cpos > 4:
		frame = frame[cpos:] + infile.read(cpos)
	else:
		stderr.write("File sync lost at byte " + str(infile.tell()) + ".\n")
		exit()
	
outfile.close()
