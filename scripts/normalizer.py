#!/usr/bin/env python2

from sys import stdout, stderr, exit
from argparse import ArgumentParser
from os.path import basename, join
from functools import partial
import logging
import csv

import numpy as np

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

def toFlt(m, xoffset=0, yoffset=0):
    ''' constructs a numpy float array from given input matrix '''
    x = len(m)
    y = len(m[xoffset])
    ary = np.full((x-xoffset, y-yoffset), float('-inf'))
    for i in xrange(len(ary)):
        for j in xrange(len(ary[i])):
            if m[i][j] != 'NULL':
                ary[i][j]  = float(m[i + xoffset][j + yoffset]) 
    return ary

def computeAvgMax(m):
    ''' compute the average of all maximum entries of all matrices '''
    matMax = float('-inf')
    diagEntrs = []

    for i in xrange(len(m)):
        for j in xrange(len(m[i])):
            if m[i][j] > matMax:
                matMax = m[i][j]
            if j == i + 1 and m[i][j] > float('-inf'):
                diagEntrs.append(m[i][j])

    # compute the average value on the diagonal above the main diagonal of
    # every chromosome
    avgAbDiag = 0
    if len(diagEntrs):
        avgAbDiag = sum(diagEntrs)/len(diagEntrs)

    # change entries to distances
    for i in xrange(len(m)):
        for j in xrange(len(m[i])):
            if m[i][j] > float('-inf'):
                m[i][j] = matMax + 1.0 - m[i][j]
            # for the main diagonal of interchromosomal maps we do not need to
            # insert artificial distances
            if i == j and len(m) == len(m[i]):
                if m[i][j] > float('-inf'):
                    m[i][j] = matMax + 1.0 - m[i][j]
                else:
                    m[i][j] = matMax + 1.0 - avgAbDiag
    return matMax


def normalizeMap(m, c):
    '''normalize a Hi-C map'''
    return c * m

def writeMap(m, orig_m, xoffset, yoffset, out):
    ''' print Hi-C map to output'''
    for k in xrange(xoffset):
        for l in xrange(len(orig_m[k])):
            out.write(orig_m[k][l])
            if l + 1 != len(orig_m[k]):
                out.write('\t')
        out.write('\n')

    for k in xrange(len(m)):
        for l in xrange(yoffset):
            out.write(orig_m[k+xoffset][l])
            out.write('\t')

        for l in xrange(len(m[0])):
            out.write(str(m[k][l]))
            # check if we have read the end of the line
            if l + 1 != len(m[0]):
                out.write('\t')
        # check if we have reached the last line
        if k + 1 != len(m):
            out.write('\n')

if __name__ == '__main__':

    parser = ArgumentParser(
            description='creates distance matrices out of Hi-C maps.')
    parser.add_argument('-o', '--out_dir', type=str, default=None, 
            help='Change output directory for normalized matrizes. By ' + \
                    'default, the directory of the input matrices is used. ')
    parser.add_argument('-x', '--row_offset', type=int, default=0, 
            help='Indicate that the <x> first rows of the matrix are ' + \
                    'row headers')
    parser.add_argument('-y', '--column_offset', type=int, default=0, 
            help='Indicate that the <y> first columns of the matrix are ' + \
                    'column names')
    parser.add_argument('map', metavar='Hi-C map', type=str, nargs='+', 
            help='path of all Hi-C maps to consider for the normalization')

    args = parser.parse_args()

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.ERROR)
    ch.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    LOG.addHandler(ch)


    LOG.info('begin loading Hi-C maps')
    mapList = [list(csv.reader(open(i), delimiter='\t'))  for i in args.map]
    arys = map(partial(toFlt, xoffset=args.row_offset,
        yoffset=args.column_offset), mapList)
    LOG.info('loading done')

    LOG.info('computing average distances')
    maxInMats = map(computeAvgMax, arys)
    avgMax = sum(maxInMats)/len(maxInMats) 
    LOG.info('done')

    for i in xrange(len(arys)):
        LOG.info('normalizing %s' %args.map[i])
        # get normalization constant 
        # NOTE TIZIAN: Taking the average here leads to smaller distances in
        # the matrices. I personally prefer this Normalization
        c = avgMax / maxInMats[i]
        m = normalizeMap(arys[i], c)
        out_file = '%s.dmat' %args.map[i]
        if args.out_dir:
            out_file = join(args.out_dir, basename(out_file))

        LOG.info('writing to %s.dmat' %out_file)
        out = open(out_file, 'w')
        writeMap(m, mapList[i], args.row_offset, args.column_offset, out)
        out.close()

    LOG.info('DONE!')
