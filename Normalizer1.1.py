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

aryList = list()
for i in range(len(mapList)):
	entries = []
	diagEntrs = []
	#a dictionary to store all values of contact counts we observe
	cCDict = {}
	#this is only to simplify notations
	mAp = mapList[i]
	ary = np.empty(mapList[i].shape)

	for i in range(len(mAp)):
		for j in range(len(mAp[i])):
			if mAp[i][j] != "NULL":
				entries.append(float(mAp[i][j]))
				if j == i + 1:
					diagEntrs.append(float(mAp[i][j]))

	#Compute the average value on the diagonal above the main diagonal of every chromosome
	avgAbDiag = 0.0

	for entr in diagEntrs:
		avgAbDiag += entr

	avgAbDiag = avgAbDiag / float(len(diagEntrs))

	matMax = max(entries)

	#change entries to distances
	for i in range(len(mAp)):
		for j in range(len(mAp[i])):
			if mAp[i][j] != "NULL":
				ary[i][j] = matMax + 1.0 - float(mAp[i][j])

				#we cast the distance to iteger to have enough identical values
				iv = int(ary[i][j])

				#fill cCDict
				if iv in cCDict:
					cCDict[iv] += 1.0
				else:
					cCDict[iv] = 1.0
			if i == j:
				ary[i][j] = matMax + 1.0 - avgAbDiag


	#get list for x and y values
	x = []

	for cC in cCDict:
		x.append(cC)

	x.sort()

	y = []

	for cCV in x:
		y.append(cCDict[cCV])

	#convert y into a numpy array
	y = np.array(y)

	est = TheilSenRegressor(random_state=42).fit(np.array(x)[:, np.newaxis], np.log(y))

	#calculate the actual coefficients
	x = [0.0, 1.0]

	y = TheilSenRegressor.predict(est, np.array(x).reshape(2, 1))

	coefList.append([(y[1] - y[0]) / (x[1] - x[0]), y[0]])

	
	aryList.append(ary)
#compute average coefficients
cZero = 0.0
cOne = 0.0

for coeffs in coefList:
	cZero += coeffs[0]
	cOne += coeffs[1]
	print 'Identified coefficients: %s' %str(coeffs)

#print coeffs
#exit()

numCoeffs = float(len(coefList))

overallCzero = cZero / numCoeffs
overallCone = cOne / numCoeffs

print "overall coefficients: zero: ", str(overallCzero), " one: ", str(overallCone)

print "Begin normalization"

#Normalize each map
for i in range(len(mapList)):
	currCzero = coefList[i][0] / overallCzero #/ coefList[i][0]
	currCone = overallCone - coefList[i][1]

	print "Coefficients we normalize with: zero: ", str(currCzero), " one: ", str(currCone)

	for k in range(len(mapList[i])):
		for l in range(len(mapList[i][k])):
			#we do not consider null entries
			if mapList[i][k][l] != "NULL":
				#calculate normalized value
#				nV = math.exp(currCzero * math.log(aryList[i][k][l]) + currCone)
				nV = currCzero * float(aryList[i][k][l]) + currCone

				aryList[i][k][l] = nV
#				#store new value in map - we do not want negative distances
#				if nV < 0.0:
#					mapList[i][k][l] = "0.0"
#				else:
#					mapList[i][k][l] = str(nV)

print "Normalisation done"

#saving the matrices
for i in range(len(mapList)):
	mapfile = open(argus.M[i] + ".dmat", "w")

	for k in range(len(mapList[i])):
		for l in range(len(mapList[i][0])):
			if mapList[i][k][l] != 'NULL':
				mapfile.write('%f' %aryList[i][k][l])
			else:
				mapfile.write('NULL')

			#check if we have read the end of the line
			if l + 1 != len(mapList[i][0]):
				mapfile.write("\t")
		#check if we have reached the last line
		if k + 1 != len(mapList[i]):
			mapfile.write("\n")
	mapfile.close()
