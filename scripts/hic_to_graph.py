#!/usr/bin/env python2

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from sys import stderr, stdout, exit, maxint
from itertools import product, combinations
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
        doSequential, out):
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

        isSymmetric = x_segments == y_segments

        if doSequential and not isSymmetric:
            LOG.warning(('Sequential mode, segments of the X and Y axis ' + \
                    'are different, thus skipping file %s') %f)
            continue

        if LOG.level == logging.DEBUG:
            inter_dist = list()
            intra_dist = list()

        last_gene = -1
        # segment counter
        i = 0
        for line in data:
            # continue with next line if no gene is associated with this
            # segment
            if not ysegs2gene[i]:
                i += 1
                continue

            if doSequential:
                std_w = 0
                if not line[i+xoffset] or line[i+xoffset] == 'NULL':
                    LOG.error(('Matrix entry at position (%s, %s) is empty' + \
                        ' or NULL, unable to assign distances to genes ' + \
                        'within segment %s') %(i, i, y_segments[i]))
                else:
                    std_w = float(line[i+xoffset])

                gi = genes[ysegs2gene[i][0]]
                if last_gene >= 0 and genes[xsegs2gene[last_gene][-1]][0] == \
                        gi[0]:

                    jp = last_gene+xoffset
                    # XXX naive heuristic for shortest path
                    # go along the x-axis, add standard distance (this is
                    # written in the diagonal) as long as off-diagonal entries
                    # are empty or NULL 
                    while jp < i+xoffset and (not line[jp] or line[jp] == \
                            'NULL'):
                        wp += std_w 
                        jp += 1

                    if jp != i + xoffset:
                        wp += float(line[jp])
                    
                    # XXX measure distance from start to start; this has been
                    # empirically studied on an arabidopsis dataset normalized
                    # with HOMER and produced a good match between intra and
                    # inter-segmental distances, whereas start to end would
                    # produce make intersegmental distances significantly
                    # shorter than intrasegmental counterparts.
                    d = abs(gi[1]-genes[xsegs2gene[last_gene][-1]][1])
                    w = wp/(y_segments[i][1] - x_segments[last_gene][1])

                    out.write(('\tedge [\n\tsource %s\n\ttarget ' + \
                            '%s\n\tweight %.4f\n\t]\n') %(xsegs2gene[
                                last_gene][-1], ysegs2gene[i][0], w*d))
                    if LOG.level == logging.DEBUG:
                        inter_dist.append(w)


                wp = std_w/(y_segments[i][2] - y_segments[i][1] + 1)
                if LOG.level == logging.DEBUG:
                    intra_dist.append(wp)
                x = ysegs2gene[i][0]
                for y in ysegs2gene[i][1:]:
                    w = wp * abs(genes[y][1] - genes[x][1])
                    out.write(('\tedge [\n\tsource %s\n\ttarget ' + \
                            '%s\n\tweight %.4f\n\t]\n') %(x, y, w))
                    x = y
                last_gene = i
            else:
                jps = xrange(xoffset, len(line))
                if isSymmetric:
                    # if symmetric, then only uses upper triangle of the
                    # distance matrix
                    jps = xrange(i+xoffset, len(line))
                for jp in jps:
                    j = jp-xoffset

                    if line[jp] and line[jp] != 'NULL' and xsegs2gene[j]:
                        wp = float(line[jp])
                        w = wp
                        if x_segments[j][0] == y_segments[i][0]:
                            # generally do comparison within the same segment
                            if i + 1 == j and LOG.level == logging.DEBUG:
                                inter_dist.append(wp/(x_segments[j][1] - \
                                        y_segments[i][1]))
                            if x_segments[i] == y_segments[j]:
                                w = delta
                            elif x_segments[j][1] > y_segments[i][1]:
                                d = abs(genes[xsegs2gene[j][-1]][1] -
                                        genes[ysegs2gene[i][0]][1])
                                # if segment distance is measured extremity to
                                # extremity we need to add +1 (same holds for
                                # gene distance in the line above)
                                #w = w/(x_segments[j][2]-y_segments[i][1]+1) * d
                                w = w/(x_segments[j][1]-y_segments[i][1]) * d
                            else:
                                d = abs(genes[ysegs2gene[i][-1]][1] -
                                        genes[xsegs2gene[j][0]][1])
                                # if segment distance is measured extremity to
                                # extremity we need to add +1 (same holds for
                                # gene distance in the line above)
                                #w = w/(x_segments[i][2]-y_segments[j][1]+1) * d
                                w = w/(y_segments[i][1]-x_segments[j][1]) * d
                        if w  > delta:
                            continue

                        p = product(xsegs2gene[j], ysegs2gene[i])
                        if x_segments[j][0] != y_segments[i][0]:
                            for x, y in p:
                                out.write(('\tedge [\n\tsource %s\n\ttarget '+\
                                        '%s\n\tweight %.4f\n\t]\n') %(x, y,
                                            wp))
                        else:
                            scale_f = x_segments[j][1] - y_segments[i][1]
                            if isSymmetric and i == j:
                                p = combinations(xsegs2gene[j], 2)
                                scale_f = x_segments[j][2] - y_segments[i][1]+1
                                if LOG.level == logging.DEBUG:
                                    intra_dist.append(wp/scale_f)

                            w = wp/scale_f

                            for x, y in p:
                                d = genes[x]!=genes[y] and \
                                        abs(genes[x][1]-genes[y][1]) or 0
                                if w*d < delta:
                                    out.write(('\tedge [\n\tsource %s\n\ttarget '+\
                                            '%s\n\tweight %.4f\n\t]\n') %(x, y,
                                                w*d))
            i += 1 
        if LOG.level == logging.DEBUG:
            LOG.debug('INTER_DISTS\t%s' %','.join(map(str, inter_dist)))
            LOG.debug('INTRA_DISTS\t%s' %','.join(map(str, intra_dist)))

                    
    # Finish the output file
    out.write(']')


if __name__ == '__main__':

    # setting up the argument parser
    parser = ArgumentParser(description='Produces a file in graphml format ' + \
            'for the gene teams algorithm', formatter_class=ADHF)
    parser.add_argument('-d', '--delta', default=maxint, type=float, 
            help='Ignore matrix values larger than DELTA to build the graph')
    parser.add_argument('-f', '--format', default = HIC_DATA_FORMATS[0], type=str,
            choices = HIC_DATA_FORMATS, 
            help='Format of the input Hi-C maps files')
    parser.add_argument('-s', '--sequential', action='store_true', 
            help='Construct graph using only distances along the diagonal, ' + \
                    'representing the linear genome sequence.')
    parser.add_argument('annotation_file', type=str, 
            help='Path to GFF3 file containing all gene annotations')
    parser.add_argument('homology_table', type=str, 
            help='Path to the file that contains names of homologous ' + \
                    'genes of a gene family. It should be a tab separated ' + \
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
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    LOG.addHandler(ch)


    genes = readAnnotations(open(args.annotation_file))
    homologies = readHomologies(open(args.homology_table))

    out = stdout
    parseHiCMapAndWriteGraph(args.hic_map, args.format, genes, homologies,
            args.delta, args.sequential, out)
    
