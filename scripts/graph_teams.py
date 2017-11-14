#!/usr/bin/env python2

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from itertools import combinations, chain, imap
from collections import deque
from random import choice
import networkx as nx

DEFAULT_DELTA = 1

#
# XXX monkey-patching new networkX API
#
if not hasattr(nx.Graph, 'nodes_iter'):
    setattr(nx.Graph, 'nodes_iter', nx.Graph.nodes)


def next_delta_connected(G, S, delta, source=None):
    """ traverses graph G by BFS: in each step, the function reports a
    (previously unreported) vertex s of S as long as there exists a path with
    distance <= delta from an already reported member to to s; requires that
    source must be part of S"""

    if source == None:
        source = choice(S)
        
    S = set(S)
    yield source
    visited = set((source, ))
    
    queue = deque((source, v, 0) for v in G.neighbors(source))
    while queue:
        u, v, d = queue.popleft()
        if v in visited:
            continue
        dv = d + G[u][v]['weight']
        if dv > delta:
            continue
        if v in S:
            dv = 0
            yield v
        visited.add(v)
        queue.extend((v, w, dv) for w in G.neighbors(v))


def smallMax(Ss, Gs, delta, source=None):
    """ identify smallest \delta-set in either G or H """
   
    BFS_handles = [next_delta_connected(Gs[i], Ss[i], delta, source and
        source[i] or None) for i in xrange(len(Gs))]

    Vs = [set() for _ in xrange(len(Gs))]
    ps = [True for _ in xrange(len(Gs))]

    while any(ps) and all(imap(lambda i: ps[i] or (len(Vs[i]) == len(Ss[i])), \
            xrange(len(Gs)))):

        for i in xrange(len(Gs)):
            try: 
                v = next(BFS_handles[i])
                Vs[i].add(v)
            except StopIteration:
                ps[i] = False

    # search for any maximal \delta-set that is smaller than its input set...
    for i in xrange(len(Gs)-1):
        if not ps[i] and len(Vs[i]) != len(Ss[i]):
            return i, Vs[i] 
    #... if none is found, return first set by default
    return len(Gs)-1, Vs[len(Gs)-1]


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


def findDeltaTeams(Gs, delta):
    """ find all graph delta-teams in G, H """
    res = list()
 
    # common classes
    classes = set([data['class'] for _, data in Gs[0].nodes(data=True)])
    for i in xrange(1, len(Gs)):
        classes.intersection_update([data['class'] for _, data in
            Gs[i].nodes(data=True)])

    # initial element of the queue is the set of common classes with graph
    # tables
    queue = [tuple(constructClassTable(G, classes) for G in Gs)]
    while queue:
        NG = queue.pop()
        Ss = [list(chain(*NG[i].values())) for i in xrange(len(Gs))]
        i, Sp = smallMax(Ss, Gs, delta)
        if len(Ss[i]) == len(Sp):
            if all(imap(lambda x: len(x) > 1, Ss)):
                res.append(NG)
        else:
            queue.append(NG)
            queue.append(division(i, Sp, Gs[i], NG))
    return res


def division(i, S, Gi, NG):
    """ parititions each table in NG, into two separate tables, overlapping
    only at vertices with shared classes. """

    NGX = [dict() for _ in NG]

    C = set()

    for v in S:
        c = Gi.node[v]['class']
        NG[i][c].remove(v)
        if not NGX[i].has_key(c):
            NGX[i][c] = set()
        NGX[i][c].add(v)
        C.add(c)

    # split class table
    for c in C:
        # copy/move entries of NGX[j] table 
        for j in xrange(len(NG)):
            if j == i:
                # remove the entry NG[i] (where S is from) if empty
                if not NG[i][c]:
                    del NG[i][c]
            else:
                if NG[i].has_key(c) and NG[i][c]:
                    NGX[j][c] = set(NG[j][c])
                else:
                    NGX[j][c] = NG[j].pop(c)
    return NGX


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('-d', '--delta', type=int, default=DEFAULT_DELTA,
            help='max-distance threshold \delta')
    parser.add_argument('-s', '--simplify', action='store_true',
            help='Instead of finding delta-teams on original graph, ' + \
                    'create a shortest-path graph with max distance ' + \
                    'delta and subsequently identify 1-teams (output' + \
                    ' is identical or normal mode)')
    parser.add_argument('graph_file', type=str, nargs='+',
            help='input graph file 1 in GML format')
    
    args = parser.parse_args()

    graphs = list()
    for f in args.graph_file:
        G = nx.read_gml(f)
        if G.is_directed():
            G = G.to_undirected()
        if args.simplify:
            G = toSPGraph(G, args.delta)
        graphs.append(G)

    if args.simplify:
        teams = findDeltaTeams(graphs, 1)
    else:
        teams = findDeltaTeams(graphs, args.delta)

    if teams:
        out = stdout
        print >> out, 'common_classes\t%s'%('\t'.join(args.graph_file))
        for NG in teams:
            print >> out, '\t'.join(chain((';'.join(map(str, \
                    sorted(NG[0].keys()))),), (';'.join(map(str, \
                    sorted(chain(*NG.values())))) for NG in NG)))
    else:
        print >> stderr, 'No %s-teams found' %args.delta

