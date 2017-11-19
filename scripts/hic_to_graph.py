#!/usr/bin/env python2

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from sys import stderr, stdout, exit, maxint
from itertools import product
import logging
import csv
import re


PAT_DIXON = re.compile('^(?:.*\|)?(\w+):(\d+)-(\d+)(?:\|.*)?$')
PAT_HOMER = re.compile('([^-]+)-(\d+)$')


LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


# default is fist value of list
HIC_DATA_FORMATS = ['DIXON', 'HOMER']


def readHomologies(data):
    ''' reading in homology table'''

    res = dict()

    # for each gene family...
    for line in csv.reader(data, delimiter='\t'):
        # store the gene family in a dictionary with the first column being the
        # family identifier
        for h in line:
            if h:
                res[h] = line[0]
    return res


def readAnnotations(data):
    ''' reading in the annotation file '''

    genes = list()

    # for each line in the annotation file...
    for ann in data:
        # split the line by its columns
        chrx, start, end, gene_id = ann.strip().split('\t')[:4]
        start, end = int(start), int(end)

        genes.append((chrx, start, end, gene_id))

    genes.sort()

    return genes


def readXSegments(data, data_type):
    ''' parse headers of Hi-C maps to determine the corresponding segments '''
    res = list()

    if data_type == 'DIXON':
        # lines that are prefixed by '#'
        line = data.next()
        while line[0].strip().startswith('#'):
            line = data.next()

        for i in xrange(1, len(line)):
            m = PAT_DIXON.match(line[i])

            if m == None:
                LOG.fatal(('Unable to determine segment from Hi-C map ' + \
                        'header \'%s\' (data_type = %s). Exiting.') %(line[i],
                            data_type))
                exit(1)
            chrx, start, end = m.groups()
            start, end = int(start), int(end)
            if chrx.lower().startswith('chr'):
                chrx = chrx[3:]
            res.append((chrx, start, end, i-1))

    elif data_type == 'HOMER':
        line = data.next()
        prev = None
        seg_len = 0
        for i in xrange(2, len(line)):
            m = PAT_HOMER.match(line[i])
            if m == None:
                LOG.fatal(('Unable to determine segment from Hi-C map ' + \
                        'header \'%s\' (data_type = %s). Exiting.') %(line[2],
                            data_type))
                exit(1)

            cur = (m.group(1), int(m.group(2)))

            if prev != None:
                if prev[0] != cur[0]:
                    res.append((prev[0], prev[1], prev[1]+seq_len-1, i-3))
                else:
                    res.append((prev[0], prev[1], cur[1]-1, i-3))
                    seq_len = cur[1]-prev[1]
            prev = cur

        if prev != None:
            res.append((prev[0], prev[1], prev[1]+seq_len-1, len(line)-3))

    return res


def readYSegments(data, data_type):
    ''' read segments from the y-axis of the matrix of the Hi-C map file'''

    res = list()

    if data_type == 'DIXON':
        isHeader = True

        c = 0
        for line in data:
            if isHeader:
                if line.startswith('#'):
                    continue
                isHeader = False
                continue
            i = line.find('\t')
            m = PAT_DIXON.match(line[:i])
            if m == None:
                LOG.fatal(('Unable to determine segment from Hi-C map ' + \
                        'header \'%s\' (data_type = %s). Exiting.') %(
                            line[:i], data_type))
                exit(1)
            chrx, start, end = m.groups()
            start, end = int(start), int(end)
            if chrx.lower().startswith('chr'):
                chrx = chrx[3:]
            res.append((chrx, start, end, c))
            c += 1

    elif data_type == 'HOMER':
        prev = None
        seg_len = 0
        # set it to -1 because we want to output the lines that have been read
        # in the iteration before
        c = -1
        isHeader = True
        for line in data:
            if isHeader:
                isHeader = False
                continue
            i = line.find('\t')
            m = PAT_HOMER.match(line[:i])
            if m == None:
                LOG.fatal(('Unable to determine segment from Hi-C map ' + \
                        'header \'%s\' (data_type = %s). Exiting.') %(line[:i],
                            data_type))
                exit(1)

            cur = (m.group(1), int(m.group(2)))

            if prev != None:
                if prev[0] != cur[0]:
                    res.append((prev[0], prev[1], prev[1]+seq_len-1, c))
                else:
                    res.append((prev[0], prev[1], cur[1]-1, c))
                    seq_len = cur[1]-prev[1]
            prev = cur

            c += 1
        if prev != None:
            res.append((prev[0], prev[1], prev[1]+seq_len-1, c))

    return res


def mapSegs2Genes(segments, genes):

    res = [list() for _ in segments]

    srt_segs = sorted(segments)
    j  = 0
    for i in xrange(len(genes)):
        chrx, start, end, _ = genes[i]
        mid = (start+end)/2
        while j < len(segments) - 1 and (segments[j][0] < chrx or \
                (chrx == segments[j][0] and \
                abs((segments[j][1] + segments[j][2])/2-mid) > \
                abs((segments[j+1][1] + segments[j+1][2])/2-mid))):
            j += 1

        if chrx == segments[j][0] and \
                min(end, segments[j][2]) - max(start, segments[j][1]) > 0:
            res[segments[j][-1]].append(i)

    return res 


def parseHiCMapAndWriteGraph(hic_map_files, data_type, genes, homologies, delta,
        out):
    # write the header
    out.write('graph [\n\tdirected 0\n')

    xoffset = 0
    if data_type == 'DIXON':
        xoffset = 1
    elif data_type == 'HOMER':
        xoffset = 2

    for i in xrange(len(genes)):
        hclass = homologies.get(genes[i][-1], genes[i][-1])
        # create node for gene
        out.write('\tnode [\n\tid %s\n\tlabel %s\n\tclass "%s"\n\t]\n' %(i, i,
            hclass))

    
    # add edges to the graph
    for f in hic_map_files:

        y_segments = readYSegments(open(f), data_type)
        ysegs2gene = mapSegs2Genes(y_segments, genes)

        data = csv.reader(open(f), delimiter='\t')
        x_segments = readXSegments(data, data_type)
        xsegs2gene = mapSegs2Genes(x_segments, genes)

        # segment counter
        i = 0
        for line in data:
            # continue with next line if no gene is associated with this segment
            if not ysegs2gene[i]:
                i += 1
                continue
            for jp in xrange(xoffset, len(line)):
                j = jp-xoffset

                if line[jp] and line[jp] != 'NULL' and xsegs2gene[j]:
                    wp = float(line[jp])
                    w = wp
                    if x_segments[j][0] == y_segments[i][0]:
                        if x_segments[j][1] > y_segments[i][1]:
                            d = abs(genes[xsegs2gene[j][-1]][2] -
                                    genes[ysegs2gene[i][0]][1])
                            w = w/(x_segments[j][2]-y_segments[i][1]) * d
                        else:
                            d = abs(genes[ysegs2gene[i][-1]][2] -
                                    genes[xsegs2gene[j][0]][1])
                            w = w/(y_segments[i][2]-x_segments[j][1]) * d
                    if w  > delta:
                        continue

                    for x, y in product(xsegs2gene[j], ysegs2gene[i]):
                        gj = genes[x]
                        gi = genes[y]

                        w = wp
                        if gj[0] == gi[0]:
                            if gj[1] > gi[1]:
                                w = w/(x_segments[j][2]  - y_segments[i][1]) * \
                                        abs(gj[2]-gi[1])
                            else:
                                w = w/(y_segments[i][2]  - x_segments[j][1]) * \
                                        abs(gi[2]-gj[1])
                        if w > delta:
                            continue

                        out.write(('\tedge [\n\tsource %s\n\ttarget ' + \
                                '%s\n\tweight %.4f\n\t]\n') %(x, y, w))
            i += 1 
                    
    # Finish the output file
    out.write(']')


if __name__ == '__main__':

    # setting up the argument parser
    parser = ArgumentParser(description='Produces a file in graphml format ' + \
            'for the gene teams algorithm', formatter_class=ADHF)
    parser.add_argument('-d', '--delta', default=maxint, type=int, 
            help='Ignore matrix values larger than DELTA to build the graph')
    parser.add_argument('-f', '--format', default = HIC_DATA_FORMATS[0], type=str,
            choices = HIC_DATA_FORMATS, 
            help='Format of the input Hi-C maps files')
    parser.add_argument('annotation_file', type=str, 
            help='Path to GFF3 file containing all gene annotations')
    parser.add_argument('homology_table', type=str, 
            help='Path to the file that contains names of homologous ' + \
                    'genes of a gene family.  Should be a tab separated ' + \
                    'file in which each line is a gene family and each ' + \
                    'column contains the names of homologous genes of a ' + \
                    'species.')
    parser.add_argument('hic_map', type=str, nargs='+', 
            help='Path to the file that contains a Hi-C map. Should be ' + \
                    'a tab separated file consisting of three colums ' + \
                    'containing the pair of bins and their count.')

    args = parser.parse_args()

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    LOG.addHandler(ch)


    genes = readAnnotations(open(args.annotation_file))
    homologies = readHomologies(open(args.homology_table))

    out = stdout
    parseHiCMapAndWriteGraph(args.hic_map, args.format, genes, homologies,
            args.delta, out)
    
