from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
mpl.style.use('seaborn-white')

import multiprocess as mp
import numpy as np
import pandas as pd
import bioframe
import cooltools
import cooler
import argparse
import os
import re

parser = argparse.ArgumentParser(description="Do compartment analysis for group of cool files")
parser.add_argument("--pca", nargs="+", required=True, help="files that have pc1 values")
parser.add_argument("--tags", nargs="+", help="Names for the cool files")
parser.add_argument("--col", type=int, help="The column that have the PC1 values")
parser.add_argument("--outdir", default=".", help="Directory to put the results in")
args = parser.parse_args()

chromsizes = bioframe.fetch_chromsizes('mm10')
chromosomes = list(chromsizes.index)

if args.tags:
	tags = args.tags
else:
	tags = []
	for pcafile in args.pca:
		filename = os.path.basename(pcafile)
		tags.append(re.sub("\.\S+", "", filename))

conditions = tags

long_names = dict(zip(tags, tags))
pal = sns.color_palette('colorblind')
colors = {tags[i] : pal[i] for i in range(len(tags))}

eigs = dict()
for i, pcafile in enumerate(args.pca):
	cond = tags[i]
	eigs[cond] = pd.read_csv(pcafile, sep="\t")


from scipy.stats import rankdata

ncols = len(tags) - 1
nrows = 4
gs = plt.GridSpec(nrows=nrows, ncols=ncols)
plt.figure(figsize=(ncols*4, nrows*3))

condx = tags[0]

for i, condy in enumerate(tags[1:]):
	plt.subplot(gs[0,i])
	lo, hi = -2 , 2
	plt.hexbin(
	    eigs[condx]['E1'],
	    eigs[condy]['E1'],
	    vmax=50,
	)
	plt.xlabel('E1 ' + long_names[condx])
	plt.ylabel('E1 ' + long_names[condy])
	plt.gca().set_aspect(1)
	plt.xlim(lo, hi)
	plt.ylim(lo, hi)
	plt.axvline(0, c='b', lw=0.5, ls='--')
	plt.axhline(0, c='b', lw=0.5, ls='--')
	plt.plot([lo, hi], [lo, hi], c='b', lw=0.5, ls='--')
	plt.colorbar(shrink=0.6)


	plt.subplot(gs[1,i])
	mask = eigs[condx]['E1'].notnull() & eigs[condy]['E1'].notnull() 
	vx = eigs[condx]['E1'].loc[mask].values
	vy = eigs[condy]['E1'].loc[mask].values
	lo, hi = 0 , len(vx)

	plt.hexbin(
	    rankdata(vx),
	    rankdata(vy),
	    vmax=20,
	)
	plt.xlabel('E1 rank ' + long_names[condx])
	plt.ylabel('E1 rank ' + long_names[condy])
	plt.gca().set_aspect(1)
	plt.xlim(lo, hi)
	plt.ylim(lo, hi)
	plt.plot([lo, hi], [lo, hi], c='b', lw=0.5, ls='--')
	plt.colorbar(shrink=0.6)

condx = tags[-1]
for i, condy in enumerate(tags[:-1]):
	plt.subplot(gs[2,i])
	lo, hi = -2 , 2
	plt.hexbin(
	    eigs[condx]['E1'],
	    eigs[condy]['E1'],
	    vmax=50,
	)
	plt.xlabel('E1 ' + long_names[condx])
	plt.ylabel('E1 ' + long_names[condy])
	plt.gca().set_aspect(1)
	plt.xlim(lo, hi)
	plt.ylim(lo, hi)
	plt.axvline(0, c='b', lw=0.5, ls='--')
	plt.axhline(0, c='b', lw=0.5, ls='--')
	plt.plot([lo, hi], [lo, hi], c='b', lw=0.5, ls='--')
	plt.colorbar(shrink=0.6)


	plt.subplot(gs[3,i])
	mask = eigs[condx]['E1'].notnull() & eigs[condy]['E1'].notnull() 
	vx = eigs[condx]['E1'].loc[mask].values
	vy = eigs[condy]['E1'].loc[mask].values
	lo, hi = 0 , len(vx)

	plt.hexbin(
	    rankdata(vx),
	    rankdata(vy),
	    vmax=20,
	)
	plt.xlabel('E1 rank ' + long_names[condx])
	plt.ylabel('E1 rank ' + long_names[condy])
	plt.gca().set_aspect(1)
	plt.xlim(lo, hi)
	plt.ylim(lo, hi)
	plt.plot([lo, hi], [lo, hi], c='b', lw=0.5, ls='--')
	plt.colorbar(shrink=0.6)

plt.savefig(f"{args.outdir}/PC1.scatter.pdf")