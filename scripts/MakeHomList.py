#!/usr/bin/env python

import sys
import csv

mfile = open(sys.argv[1], "r")
homList = [(cols[0], cols[3]) for cols in csv.reader(mfile, delimiter='\t')]
mfile.close()

hfile = open(sys.argv[2], "r")
humanGeneDict = {}

for cols in csv.reader(hfile, delimiter='\t'):
    if cols[0] not in humanGeneDict:
        humanGeneDict[cols[0]] = cols[3]

hfile.close()

hFile = open(sys.argv[3], "w")
for ent1, ent2 in homList:
    hFile.write(humanGeneDict[ent1] + '\t' + ent2+ "\n")
