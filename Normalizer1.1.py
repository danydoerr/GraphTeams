#!/usr/bin/env python

#In this version of the normalizer the normalization is done using Theil-Sen Regression

import sys
import numpy as np
import argparse as args
from sklearn.linear_model import TheilSenRegressor
import math

parser = args.ArgumentParser(description="Creates distance matrices out of Hi-C maps.")
parser.add_argument('M', metavar='Hi-C_Map', type=str, nargs='+', help="Path of all Hi-C maps to consider for the normalization")

argus = parser.parse_args()

mapList = [np.genfromtxt(i, dtype=str) for i in argus.M]

print "Loading done"

#convert contact counts in the matrices to distances and count the number of each different contact count in the matrices
#list for the normalization coefficients of each map
coefList = []

for i in range(len(mapList)):
	entries = []
	#a dictionary to store all values of contact counts we observe
	cCDict = {}
	#this is only to simplify notations
	mAp = mapList[i]

	for i in range(len(mAp)):
		for j in range(len(mAp[i])):
			if mAp[i][j] != "NULL":
				entries.append(float(mAp[i][j]))

	matMax = max(entries)

	#change entries to distances
	for i in range(len(mAp)):
		for j in range(len(mAp[i])):
			if mAp[i][j] != "NULL":
				mAp[i][j] = str(matMax + 1.0 - float(mAp[i][j]))

				#we cast the distance to iteger to have enough identical values
				iv = str(int(float(mAp[i][j])))

				#fill cCDict
				if iv in cCDict:
					cCDict[iv] += 1.0
				else:
					cCDict[iv] = 1.0

	#get list for x and y values
	x = []

	for cC in cCDict:
		x.append(float(int(cC)))

	x.sort()

	y = []

	for cCV in x:
		y.append(cCDict[str(int(cCV))])

	#convert y into a numpy array
	y = np.array(y)

	est = TheilSenRegressor(random_state=42).fit(np.array(x)[:, np.newaxis], np.log(y))

	#calculate the actual coefficients
	x = [0.0, 1.0]

	y = TheilSenRegressor.predict(est, np.array(x).reshape(2, 1))

	coefList.append([(y[1] - y[0]) / (x[1] - x[0]), y[0]])

#compute average coefficients
cZero = 0.0
cOne = 0.0

for coeffs in coefList:
	cZero += coeffs[0]
	cOne += coeffs[1]

numCoeffs = float(len(coefList))

overallCzero = cZero / numCoeffs
overallCone = cOne / numCoeffs

print "Begin normalization"

#Normalize each map
for i in range(len(mapList)):
	currCzero = overallCzero / coefList[i][0]
	currCone = overallCone / coefList[i][1]

	for k in range(len(mapList[i])):
		for l in range(len(mapList[i][k])):
			#we do not consider null entries
			if mapList[i][k][l] != "NULL":
				#calculate normalized value
				nV = currCzero * math.log(float(mapList[i][k][l])) + currCone

				#store new value in map - we do not want negative distances
				if nV < 0.0:
					mapList[i][k][l] = "0.0"
				else:
					mapList[i][k][l] = str(nV)

print "Normalisation done"

#saving the matrices
for i in range(len(mapList)):
	mapfile = open(argus.M[i] + ".dmat", "w")

	for k in range(len(mapList[i])):
		for l in range(len(mapList[i][0])):
			mapfile.write(mapList[i][k][l])

			#check if we have read the end of the line
			if l + 1 != len(mapList[i][0]):
				mapfile.write("\t")
		#check if we have reached the last line
		if k + 1 != len(mapList[i]):
			mapfile.write("\n")
	mapfile.close()
