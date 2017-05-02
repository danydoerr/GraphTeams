#!/usr/bin/env python

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from itertools import combinations, chain
from collections import deque
import networkx as nx

DEFAULT_DELTA = 1
DEFAULT_OUT_PREFIX = 'teams'

#
# XXX monkey-patching new networkX API
#
if not hasattr(nx.Graph, 'nodes_iter'):
    setattr(nx.Graph, 'nodes_iter', nx.Graph.nodes)


def next_delta_connected(G, S, source, delta):
    """ traverses graph G by BFS, reporting vertices of S as long as any two of
    its members (which must include the source) is closer than delta """
 
    S = set(S)
    if G.node[source]['class'] not in S:
        raise Exception('Source not part of set S')
    yield source, G.node[source]['class']

    visited = set((source, ))
    
    queue = deque((source, v, 0) for v in G.neighbors(source))
    while queue:
        u, v, d = queue.popleft()
        if v in visited:
            continue
        dv = d + G[u][v]['weight']
        if dv > delta:
            continue
        s = G.node[v]['class']
        if s in S:
            dv = 0
            yield v, s
        visited.add(v)
        queue.extend((v, w, dv) for w in G.neighbors(v))


def smallMax(S, G, H, sG, sH, delta):
    """ identify smallest \delta-set in either G or H """
    
    BFS_G = next_delta_connected(G, S, sG, delta)
    BFS_H = next_delta_connected(H, S, sH, delta)
    V_G = set()
    V_H = set()
    S_G = set()
    S_H = set()
    while True:
        try: 
            v, s = next(BFS_G)
            S_G.add(s)
            V_G.add(s)
        except StopIteration:
            return 1, V_G, S_G
        try:
            v, s = next(BFS_H)
            S_H.add(s)
            V_H.add(v)
        except StopIteration:
            return 2, V_H, S_H


def constructClassTable(G, classes):
    """ constructs dictionary representing table with classes and their
    associated vertices """
    res = dict((c, set()) for c in classes)
    for v, data in G.nodes(data=True):
        c = data['class']
        if c in classes:
            res[c].add(v)
    return res


def toSPGraph(G, delta):
    """ convert to shortest-path graph """
    
    SP = nx.shortest_path_length(G, weight='weight')
    spG = nx.Graph()
    spG.add_nodes_from(G.nodes(data=True))
    
    for u in SP:
        for v, d in SP[u].items():
            if u != v and d <= delta:
                spG.add_edge(u, v, weight=1)
    return spG


def findDeltaTeams(G, H, delta):
    """ find all graph delta-teams in G, H """
    res = list()
 
    # common classes
    classes = set([data['class'] for _, data in G.nodes(data=True)])
    classes.intersection_update([data['class'] for _, data in
        H.nodes(data=True)])

    # initial element of the queue is the set of common classes with graph
    # tables
    queue = [(constructClassTable(G, classes), constructClassTable(H, classes))]
    while queue:
        NG, NH = queue.pop()
        source = next(iter(NG))
        i, V, Sp = smallMax(NG.keys(), G, H, next(iter(NG[source])),
                next(iter(NH[source])), delta)
        if len(Sp) == len(NG):
            if len(Sp) > 1:
                res.append((NG, NH))
        else:
            queue.append((NG, NH))
            queue.append(division(i, V, Sp, NG, NH))
    return res


def division(i, V, Sp, NG, NH):
    """ parititions tables NG, NH into four separate graphs and
    tables, overlapping only at vertices with shared classes. """

    # split first the table where V is from
    if i == 2:
        NH, NG = NG, NH
    # split class table
    NGX = dict()
    NHX = dict()
    for s in Sp:
        # move vertices not in V of table NG to NGX
        NGX[s] = set()
        for v in list(NG[s]):
            if v in V:
                NG[s].remove(v)
                NGX[s].add(v)
        # copy/move entries of NHX table 
        if len(NG[s]) and len(NGX[s]):
            NHX[s] = set(NH[s])
        elif len(NGX[s]):
            del NG[s]
            NHX[s] = NH.pop(s)
        else:
            del NGX[s]
    # switch again to original order
    if i == 2:
        NHX, NGX = NGX, NHX
    return NGX, NHX


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('-o', '--out_file_prefix', type=str,
            default=DEFAULT_OUT_PREFIX,
            help='Specify prefix for output files. The ending will ' + \
                    '*always* be ".1.gml" and ".2.gml", respectively')
    parser.add_argument('-d', '--delta', type=int, default=DEFAULT_DELTA,
            help='max-distance threshold \delta')
    parser.add_argument('-s', '--simplify', action='store_true',
            help='Instead of finding delta-teams on original graph, ' + \
                    'create a shortest-path graph with max distance ' + \
                    'delta and subsequently identify 1-teams (output' + \
                    ' is identical or normal mode)')
    parser.add_argument('graph_file1', type=str, 
            help='input graph file 1 in GML format')
    parser.add_argument('graph_file2', type=str, 
            help='input graph file 2 in GML format')
    
    args = parser.parse_args()
    
    G = nx.read_gml(args.graph_file1)
    H = nx.read_gml(args.graph_file2)

    if G.is_directed():
        G = G.to_undirected()
    if H.is_directed():
        H = H.to_undirected()

    if args.simplify:
        G = toSPGraph(G, args.delta)
        H = toSPGraph(H, args.delta)
        teams = findDeltaTeams(G, H, 1)
    else:
        teams = findDeltaTeams(G, H, args.delta)

    if teams:
        print >> stdout, 'common_classes\t%s\t%s'%(args.graph_file1,
                args.graph_file2)
        for NG, NH in teams:
            print '\t'.join((';'.join(map(str, NG.keys())), ';'.join(map(str,
                chain(*NG.values()))), ';'.join(map(str, chain(*NH.values())))))
    else:
        print >> stderr, 'No %s-teams found' %args.delta

