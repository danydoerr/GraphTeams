#!/usr/bin/env python

import sys

efile = open(sys.argv[1], "r")
oFile = open(sys.argv[2], "w")

#If we do not want a confidence filter, we comment this in
#for line in efile:
#	cols = line.split(',')
#
#	oFile.write(cols[5] + '\t' + cols[6] + '\t' + cols[7] + '\t' + cols[3] + "\n")

#If we want a confidence filter, we comment this in
for line in efile:
	line = line.replace("\n", "")

	cols = line.split(',')

	if cols[9] == '1':
		oFile.write(cols[5] + '\t' + cols[6] + '\t' + cols[7] + '\t' + cols[3] + "\n")
