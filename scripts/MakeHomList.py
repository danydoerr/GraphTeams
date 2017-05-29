#!/usr/bin/env python

import sys

mfile = open(sys.argv[2], "r")

homList = []

for line in mfile:
	cols = line.split(',')

	homList.append((cols[0], cols[3]))

mfile.close()

hfile = open(sys.argv[1], "r")

humanGeneDict = {}

for line in hfile:
	cols = line.split(',')

	if not cols[2] in humanGeneDict:
		humanGeneDict[cols[2]] = cols[3]

hfile.close()

hFile = open(sys.argv[3], "w")

for ent in homList:
	hFile.write(humanGeneDict[ent[0]] + '\t' + ent[1]+ "\n")
