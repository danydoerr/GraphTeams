#!/usr/bin/env python2

from sys import stdout, stderr, exit
from argparse import ArgumentParser
from os.path import basename, join
from functools import partial
from math import isinf
import logging
import csv

import numpy as np

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

NORMALIZATION_TYPES = ['MAX', 'AVG']



def toFlt(m, row_offset=0, column_offset=0):
    ''' constructs a numpy float array from given input matrix '''
    x = len(m)
    y = len(m[row_offset])
    ary = np.full((x-row_offset, y-column_offset), float('-inf'))
    for i in xrange(len(ary)):
        for j in xrange(len(ary[i])):
            if m[i][j] != 'NULL':
                ary[i][j]  = float(m[i + row_offset][j + column_offset]) 
    return ary


def toDist(m):
    ''' compute the average of all maximum entries of all matrices '''

    new_m = m - np.min(m[m > float('-inf')])
    matMax = np.max(new_m)
    diagEntrs = []

    for i in xrange(len(new_m)-1):
        if new_m[i][i+1] > float('-inf'):
            diagEntrs.append(new_m[i][i+1])

    # compute the average value on the diagonal above the main diagonal of
    # every chromosome
    avgAbDiag = 0
    if len(diagEntrs):
        avgAbDiag = sum(diagEntrs)/len(diagEntrs)

    new_m = matMax - new_m

    # increase resolution of distances for the main diagonal of
    # intrachromosomal maps 
    if len(new_m) == len(new_m[0]):
        # change entries to distances
        for i in xrange(len(new_m)):
            if isinf(new_m[i][i]):
                new_m[i][i] = matMax - avgAbDiag
    return new_m


def normalizeMap(m, c):
    ''' normalize a Hi-C map '''
    return c * m


def writeMap(m, orig_m, row_offset, column_offset, out):
    ''' print Hi-C map to output '''
    for k in xrange(row_offset):
        for l in xrange(len(orig_m[k])):
            out.write(orig_m[k][l])
            if l + 1 != len(orig_m[k]):
                out.write('\t')
        out.write('\n')

    for k in xrange(len(m)):

        for l in xrange(column_offset):
            out.write(orig_m[k+row_offset][l])
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
    parser.add_argument('-x', '--column_offset', type=int, default=0, 
            help='Indicate that the <x> first columns of the matrix are ' + \
                    'column names')
    parser.add_argument('-y', '--row_offset', type=int, default=0, 
            help='Indicate that the <y> first rows of the matrix are ' + \
                    'row headers')
    parser.add_argument('-t', '--type', type=str, choices=NORMALIZATION_TYPES,
            default=NORMALIZATION_TYPES[0],
            help='Choose between different types of normalization')
    parser.add_argument('map', metavar='Hi-C map', type=str, nargs='+', 
            help='path of all Hi-C maps to consider for the normalization')

    args = parser.parse_args()

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    LOG.addHandler(ch)


    LOG.info('begin loading Hi-C maps..')
    mapList = [list(csv.reader(open(i), delimiter='\t'))  for i in args.map]
    arys = map(partial(toFlt, row_offset=args.row_offset,
        column_offset=args.column_offset), mapList)
    LOG.info('...done')

    LOG.info('transforming into distance matrix...')
    arys_dist = map(toDist, arys)
    LOG.info('...done')

    c = [0] * len(arys_dist)
    if args.type == 'MAX':
        maxInMats = map(np.max, arys_dist)
        avgMax = sum(maxInMats)/len(maxInMats)
        c = map(lambda x: avgMax/x, maxInMats)
    elif args.type == 'AVG':
        meanInMats = map(lambda x: np.mean(x[x < float('inf')]), arys_dist)
        avgMean = sum(meanInMats)/len(meanInMats)
        c = map(lambda x: avgMean/x, meanInMats)
    else:
        raise Exception, 'Normalization type %s implemented' %args.type

    for i in xrange(len(arys_dist)):
        LOG.info('normalizing %s' %args.map[i])
        # get normalization constant 
        # NOTE TIZIAN: Taking the average here leads to smaller distances in
        # the matrices. I personally prefer this Normalization
        m = normalizeMap(arys_dist[i], c[i])
        LOG.info('\tfinal average distance\t%s' %np.mean(m))
        LOG.info('\tfinal maximum distance\t%s' %np.max(m))

        out_file = '%s.dmat' %args.map[i]
        if args.out_dir:
            out_file = join(args.out_dir, basename(out_file))

        LOG.info('writing to %s.dmat' %out_file)
        out = open(out_file, 'w')
        writeMap(m, mapList[i], args.row_offset, args.column_offset, out)
        out.close()

    LOG.info('DONE!')
