#!/usr/bin/python

'''
Created on Jan 17, 2010

@author: wibble
'''
from struct import unpack,pack
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

if len(argv) < 2 or '-h' in argv or '--help' in argv or '/h' in argv:
    help()
    exit()

stderr.write('Input: Binary from ' + argv[1] + '.\n')
infile = open(argv[1],'r')

if len(argv) < 3:
    outfile = stdout
    stderr.write('Output: ASCII to standard output.\n')
else:
    outfile = open(argv[2],'w')
    stderr.write('Output: Binary to ' + argv[2] + '.\n')
    tfo = open(argv[2]+'.time','w')
    stderr.write("Time: ASCII to " + argv[2] + ".time.\n")

ptpfifo = infile.read(10000)
tmfifo = ""
dmfifo = ""

ptpfr_num = 0
colfr_num = 0
majfr_num = 0

while len(ptpfifo) >= 814:
	
	# sync to PTP
	start = ptpfifo.find("\x30\x03\x03\x30\x00\x01\x26\x03")
	end = ptpfifo.find("\x30\x03\x03\x30\x00\x01\x26\x03", start+1)
	
	ptp_frame_len = end-start
	
	if start != 0:
		print "Discontinuity?"
	
	# should always be a full PTP frame, extract
	if ptp_frame_len != 814:
		print "PTP Frame #{0}, bad size {1}".format(ptpfr_num, ptp_frame_len)
		
		ptpfifo = ptpfifo[end:]
		ptpfifo += infile.read(ptp_frame_len)
		
		continue
		
	ptpfr_num += 1

	# save latest timestamp
	tstamp = unpack('>Q',"\x00\x00"+ptpfifo[8:14])[0]
	sec = ((tstamp>>10) & 0x7FFFFFF)/ 1000.0 # 27-bit milliseconds-of-day
	sec += (tstamp & 0x3FF)*0.000001 # 10-bit microseconds-of-millisecond
	sec -= 35350.234 # Subtract launch time to get seconds-since-launch
	
	# store to TM fifo, trim off extracted section
	tmfifo += ptpfifo[start:end]	
	ptpfifo = ptpfifo[end:]
	ptpfifo += infile.read(end)

	# sync to TM (major frame)
	start = tmfifo.find("\xFE\x6B\x28\x40")
	end = tmfifo.find("\xFE\x6B\x28\x40",start+1)
		
#	if start != 0:
#		print "Major frame sync not at PTP[0], frame #{0}.".format(tm_frame_num)
	
	tm_frame_len = end-start
	
	if end > 0: # we have two TM syncs
		if tm_frame_len >= 800: # we have a full TM frame
			dmfifo += tmfifo[start+6:start+394] + tmfifo[start+398:start+800]
			
			majfr_num += 1
	
		else:
			print("TM Frame {0} @ {1} too short ({2}), discarding...".format(majfr_num, sec, tm_frame_len))
			
		tmfifo = tmfifo[end:]
	
	# sync to Colonel Frame
	start = dmfifo.find("Dartmouth College")
	end = dmfifo.find("Dartmouth College", start+1)
	
	dm_frame_len = end-start
	
	if end > 0: # we have two CF syncs
		if dm_frame_len >= 131070:
			# we have at least a full Colonel Frame (2x65535 bytes)
		
			wstart = start + 50
			wend = end if dm_frame_len == 131070 else start+131070
		
			outfile.write(dmfifo[wstart:wend])
			tfo.write("{0}\n".format(sec))
			
			colfr_num += 1
			
		else:
			print("Colonel Frame {0} @ {1} too short ({2}), discarding...".format(colfr_num, sec, dm_frame_len))
			
		dmfifo = dmfifo[end:]
			
					
print("{0} ptp, {1} maj, {2} col.\n".format(ptpfr_num, majfr_num, colfr_num))	

outfile.close()
tfo.close()