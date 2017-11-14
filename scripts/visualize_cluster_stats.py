#!/usr/bin/env python2

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from matplotlib import pylab as plt

DEFAULT_DELTA = 1
COLOR1 = [0.251,0.498,0.498]
COLOR2 = [0.616,0.306,0.537]
COLOR3 = [0.694,0.792,0.396]
COLOR4 = [0.831,0.604,0.416]

def visualize_stats(data, out):
    f, axs = plt.subplots(2, 2, sharey=False, sharex=False)
    
    X = plt.arange(data.shape[0]) 
    width = 0.4
    title_fontsize = 10
    label_fontsize = 6

    # number of detected clusters
    bars3dsize = axs[0][0].bar(X, data[:, 1], width, color=COLOR1)
    bars1dsize = axs[0][0].bar(X+width, data[:, 2], width, color=COLOR2)
    axs[0][0].legend((bars3dsize[0], bars1dsize[0]), ('3D gene clusters', \
            '1D gene clusters'), fontsize=title_fontsize)
    axs[0][0].set_xticks(X+width/2)
    axs[0][0].set_xticklabels(map(int, data[:, 0]))
    axs[0][0].set_ylabel('count')
#    axs[0][0].set_xlabel(r'$\delta$')
    axs[0][0].set_title('number of discovered clusters', fontsize=title_fontsize)

    # average size of detected clusters
    bars3dsize = axs[0][1].bar(X, data[:, 3], width, color=COLOR1)
    bars1dsize = axs[0][1].bar(X+width, data[:, 4], width, color=COLOR2)
    axs[0][1].legend((bars3dsize[0], bars1dsize[0]), ('3D gene clusters', \
            '1D gene clusters'), fontsize=title_fontsize)
    axs[0][1].set_xticks(X+width/2)
    axs[0][1].set_xticklabels(map(int, data[:, 0]))
    axs[0][1].set_ylabel('avg. size')
#    axs[0][1].set_xlabel(r'$\delta$')
    axs[0][1].set_title('average size of discovered clusters',
            fontsize=title_fontsize)
    axs[0][1].set_yscale('log')

    # runtime of GraphTeams
    axs[1][1].scatter(X+width/2, data[:, 5], color=COLOR1)
#    axs[1][0].legend(('1D gene clusters', ), fontsize=title_fontsize)
    axs[1][1].set_xticks(X+width/2)
    axs[1][1].set_xticklabels(map(int, data[:, 0]))
    axs[1][1].set_ylim([0, max(data[:,5])+10])
    axs[1][1].set_ylabel('time in mins')
    axs[1][1].set_xlabel(r'$\delta$')
    axs[1][1].set_title('computation time of 3D clusters', fontsize=title_fontsize)

    # spatial gain 3D clusters
    axs[1][0].bar(X+(width/2), data[:, 6], width, color=COLOR1)
#    axs[1][1].legend(('3D gene clusters', ), fontsize=title_fontsize)
    axs[1][0].set_xticks(X+width/2)
    axs[1][0].set_xticklabels(map(int, data[:, 0]))
    axs[1][0].set_ylabel('avg. number of gained genes')
    axs[1][0].set_xlabel(r'$\delta$')
    axs[1][0].set_title('spatial gain of 3D clusters', fontsize=title_fontsize)
    axs[1][0].set_yscale('log')

    for axsx in axs:
        for axsy in axsx:
            for label in axsy.get_xmajorticklabels():
                label.set_rotation(90)
                label.set_fontsize(label_fontsize)
                #label.set_horizontalalignment('right')
            for label in axsy.get_ymajorticklabels():
                label.set_fontsize(label_fontsize)

    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    f.savefig(out, format='pdf')

if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('cluster_statistics', type=str, 
            help='TAB-separated file containing statistical assessment ' + \
                    'from 3D gene cluster evaluation')
    args = parser.parse_args()
   
    out = stdout
    data = plt.loadtxt(args.cluster_statistics)
    visualize_stats(data, out)
