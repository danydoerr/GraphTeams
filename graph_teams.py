#!/usr/bin/env python

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from itertools import combinations, chain
import networkx as nx

DEFAULT_DELTA = 1
DEFAULT_OUT_PREFIX = 'teams'

def pruneToClasses(G, classes, delta=-1):
    """ remove all nodes from graph that do not belong to given set of classes;
    in doing so, establish direct edges between neighbors of removed vertices
    """
    for v in G.nodes():
        if v not in classes:
            for u, w in combinations(G.neighbors(v), 2):
                if not G.has_edge(u, w) and (delta < 0 or
                        G[u][v]['weight']+G[v][w]['weight'] <= delta):
                    G.add_edge(u, w, weight=G[u][v]['weight']+G[v][w]['weight'])
    return G

def constructClassTable(G):
    """ constructs dictionary representing table with classes and their
    associated vertices """
    res = dict()
    for v, data in G.nodes(data=True):
        c = data['class']
        if not res.has_key(c):
            res[c] = set()
        res[c].add(v)
    return res


def findTeams(G, H, delta):
    """ find all graph teams in G, H with max-distance delta between any two
    connected nodes """
    res = list()
 
    # initial element of the queue are both graphs with their class tables
    queue = [(G, constructClassTable(G), H, constructClassTable(H))]
    while queue:
        G, NG, H, NH = queue.pop()
        # do not output graphs containing single vertices
        if nx.is_connected(G) and nx.is_connected(H):
            if G.size() or H.size():
                res.append((G, H))
        else:
            GX, NGX, HX, NHX = division(G, NG, H, NH, delta)
            queue.extend(((GX, NGX, HX, NHX), (G, NG, H, NH)))
    return res


def division(G, NG, H, NH, delta):
    """ parititions graphs G, H, and talbes NG, NH into four separate graphs and
    tables, overlapping only at vertices with shared classes. NOTE: function
    will not terminate if both G and H are connected! """

    if nx.is_connected(G):
        HX, NHX, GX, NGX = division(H, NH, G, NG, delta)
    else:
        # get an arbitrary connected component from G
        # XXX workaround for networkx 2.0 change of API
        v = None
        if hasattr(G, 'nodes_iter'):
            v = next(G.nodes_iter())
        else:
            # nodes() returns now a generator in NetworkX version 2.0
            v = next(G.nodes())
        C = nx.node_connected_component(G, v)
        GX = G.subgraph(C)
        G.remove_nodes_from(C)
        # split class table
        NGX = dict()
        NHX = dict()
        # extract NGX table from NG
        for v, data in GX.nodes(data=True):
            f = data['class']
            if NG.has_key(f):
                if not NGX.has_key(f):
                    NGX[f] = set()
                NG[f].remove(v)
                NGX[f].add(v)
                # remove column if no longer used
                if not NG[f]:
                    del NG[f]
        # copy/move entries of NHX table 
        for f in NGX.keys():
            if not NG.has_key(f):
                NHX[f] = NH.pop(f)
            else:
                NHX[f] = set(NH[f])

        HX = H.subgraph(chain(*NHX.values()))
        # capture distances between vertices with paths not in HX
        boundary = set()
        for (u, v) in nx.edge_boundary(H, HX.nodes()):
            if HX.has_node(u):
                boundary.add(u)
            else:
                boundary.add(v)
        for u, v in combinations(boundary, 2):
            if not HX.has_edge(u, v):
                w = nx.dijkstra_path_length(H, u, v, weight='weight')
                if w <= delta:
                    HX.add_edge(u, v, weight=w)
        # remove HX from H
        H.remove_nodes_from(chain(*(NHX[f] for f in NHX.keys() if not
            NH.has_key(f))))
    return GX, NGX, HX, NHX

if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('-o', '--out_file_prefix', type=str,
            default=DEFAULT_OUT_PREFIX,
            help='Specify prefix for output files. The ending will ' + \
                    '*always* be ".1.gml" and ".2.gml", respectively')
    parser.add_argument('-d', '--delta', type=int, default=DEFAULT_DELTA,
            help='max-distance threshold \delta')
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

    # common classes
    classes = set([data['class'] for _, data in G.nodes(data=True)])
    classes.intersection_update([data['class'] for _, data in
        H.nodes(data=True)])
    pruneToClasses(G, classes) 
    pruneToClasses(H, classes) 

    # remove edges that don't satisfy delta
    for u, v, data in G.edges(data=True):
        if data['weight'] > args.delta:
            G.remove_edge(u, v)
    for u, v, data in H.edges(data=True):
        if data['weight'] > args.delta:
            H.remove_edge(u, v)

    teams = findTeams(G, H, args.delta)
    out1 = open('%s.1.gml' %args.out_file_prefix, 'w')
    out2 = open('%s.2.gml' %args.out_file_prefix, 'w')

    for GX, HX in teams:
        nx.write_gml(GX, out1)
        nx.write_gml(HX, out2)


