import multiprocess as mp
import numpy as np
import pandas as pd
import bbi
import argparse
import os
import re

from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
mpl.style.use('seaborn-white')


def plotHeatmap(clusters, cnames, bw, conditions, outdir, swap=None, nbins=50, flank=500000, vmin=-1, vmax=1, prefix='PC1'):
	stacks = {}
	for j, region in enumerate(cnames):
		peaks = pd.read_csv(clusters[j], names=['chrom', 'start', 'end'], sep="\t")
		for i, cond in enumerate(conditions):
			mids = (peaks['start'] + peaks['end']) // 2
			bwfile = bw[i]
			stacks[cond] = bbi.stackup(bwfile, peaks['chrom'], mids - flank, mids + flank, bins=nbins)

			if swap:
				for i in range(len(stacks[cond])):
					pcarray = stacks[cond][i]
					if np.nanmean(pcarray[0:25]) < np.nanmean(pcarray[25:50]):
						stacks[cond][i] = np.flip(pcarray)

		gs = GridSpec(nrows=3, ncols=len(conditions), height_ratios=[15, 2, 0.5], hspace=0)
		plt.figure(figsize=(3 * len(conditions), 10))

		X = stacks[conditions[0]]
		idx = np.argsort(X[:, X.shape[1] // 2])
		x = np.linspace(-flank / 1e6, flank / 1e6, nbins)
		cmap = plt.cm.get_cmap('coolwarm')
		cmap.set_bad('#777777')
		im_opts = dict(vmin=vmin, vmax=vmax, extent=[-flank / 1e6, flank / 1e6, len(peaks), 0], cmap=cmap)

		for i, name in enumerate(stacks):
			# heatmap
			ax = ax1 = plt.subplot(gs[0, i])
			X = stacks[name]
			img = ax.matshow(X[idx, :], **im_opts, rasterized=True)
			ax.axvline(0, c='grey', lw=0.5)
			ax.grid('off')
			ax.set_aspect('auto')
			ax.set_title(name)
			if i > 0:
				ax.yaxis.set_visible(False)

			# summary
			ax = plt.subplot(gs[1, i], sharex=ax1)
			ax.axhline(0, c='#777777', lw=1, ls='--')
			ax.plot(x, np.nanmean(stacks[name], axis=0), c='k', lw=2)
			ax.set_xlim(-flank / 1e6, flank / 1e6)
			ax.xaxis.set_visible(False)
			ax.set_ylim(vmin, vmax)
			if i > 0:
				ax.yaxis.set_visible(False)

			# color bar
			cax = plt.subplot(gs[2, i])
			cb = plt.colorbar(img, cax=cax, orientation='horizontal')
			cb.locator = mpl.ticker.MaxNLocator(nbins=3)
			cb.update_ticks()
		plt.savefig(f'{outdir}/{region}.{prefix}.pdf')


def main():
	parser = argparse.ArgumentParser(description="Plot target region PC1 values")
	parser.add_argument("-bed", nargs="+", required=True, help="bed files of target region")
	parser.add_argument("-names", nargs="+", help="names of the bed regions")
	parser.add_argument("-bw", nargs="+", required=True, help="bigwig files of PC1")
	parser.add_argument("-conds", nargs="+", help="conditions of the bw files")
	parser.add_argument("-outdir", default=".", help="the direcotry to put the output files")
	parser.add_argument("-swap", default=False, action='store_true', help="swap the region so that left is larger than right")
	parser.add_argument("-flank", default=500000, type=int, help="The length of flanking region")
	parser.add_argument("-bins", default=50, type=int, help="Number of bins of the target region")
	parser.add_argument("-vmin", default=-1, type=float, help="The min value of heatmap")
	parser.add_argument("-vmax", default=1, type=float, help="The max value of heatmap")
	parser.add_argument("-prefix", default="PC1", help="Prefix for the heatmap")
	args = parser.parse_args()

	if args.names:
		names = args.names
	else:
		names = []
		for bedfile in args.bed:
			filename = os.path.basename(bedfile)
			names.append(re.sub("\.\S+", "", filename))

	conds = []
	if args.conds:
		conds = args.conds
	else:
		conds = []
		for bwfile in args.bw:
			filename = os.path.basename(bwfile)
			conds.append(re.sub("\.\S+", "", filename))

	plotHeatmap(args.bed, names, args.bw, conds, args.outdir, args.swap, nbins=args.bins, flank=args.flank, vmin=args.vmin, vmax=args.vmax, prefix=args.prefix)


if __name__ == '__main__':
	main()
