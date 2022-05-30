#!/usr/bin/env python

import os
import sys
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import argparse

def tplot(rdf:pd.DataFrame, windowsize:int):
    """creates plot of contigs in order of size from largest to smallest with mapped telomeric read depths"""
    # list of contigs in order length
    tigs = list(rdf[rdf['pos'] > 40].groupby(['tig'], sort=False)['pos'].max().sort_values(ascending=False).keys())
    fig,ax = plt.subplots(len(tigs),1,figsize=(8,20),sharey = True, sharex = True)
    #maxx = rdf['pos'].max() * windowsize
    maxy = rdf['depth'].max()
    
    for i,tig in enumerate(tigs):
        tdf = rdf[rdf['tig'] == tig]
        ax[i].bar(x=tdf['pos'] * windowsize / 1000,height=tdf['depth'] ,width= windowsize / 500,color = 'orange')
        ax[i].plot((0,tdf['pos'].max() * windowsize / 1000),(0,0),color = 'blue',linewidth = 4)
        if not i == len(tigs) -1:
            ax[i].axis('off')
        else:
            ax[i].spines["top"].set_visible(False)
            ax[i].spines["right"].set_visible(False)
            ax[i].spines["left"].set_visible(False)
            ax[i].get_yaxis().set_visible(False)
            ax[i].set_xlabel('position (Mbp)', fontsize = 14)
    xticklabels = [ int(k) for k in ax[-1].get_xticks() ]
    ax[-1].set_xticklabels(xticklabels, fontsize = 14)
    return fig,ax

def main():
    # Define command line arguments
    parser = argparse.ArgumentParser("""This function take a depth file from samtools depth and plots depths within windows""")
    parser.add_argument('--depthfile', help = 'output file from "samtools depth"', required = True)
    parser.add_argument('--windowsize', help = 'define the window size to bin read counts', default = 10000, type = int)
    parser.add_argument('--output', help = 'file name to write plot, file type determined by suffix (e.g., species1telo.png)', required = True)
    args = parser.parse_args()

    # reading depthfile from samtools and reads it into a dictionary 
    readsdict = {}
    for line in open(args.depthfile):
        fields = line.strip().split('\t')
        pos = int(int(fields[1])/args.windowsize)
        seq = fields[0]
        dictkey = seq + ':' + str(pos)
        depth = int(fields[2])
        if not dictkey in readsdict:
            readsdict[dictkey] = [seq,pos,0]
        readsdict[dictkey][2] += depth

    # parses dictionary it into a dataframe
    rdf = pd.DataFrame().from_dict(readsdict,orient='index')
    rdf.columns = ['tig','pos','depth']
    rdf.head()

    fig,ax = tplot(rdf,args.windowsize)
    fig.savefig(args.output)
    
if __name__ == '__main__':
    main()