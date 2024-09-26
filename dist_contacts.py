# import core packages
import warnings
warnings.filterwarnings("ignore")
from itertools import combinations

# import semi-core packages
import matplotlib.pyplot as plt
from matplotlib import colors
plt.style.use('seaborn-poster')
import seaborn as sns
import numpy as np
import pandas as pd

# import open2c libraries
import bioframe

import cooler
import cooltools
import cooltools.expected

import os
import re
import argparse

parser = argparse.ArgumentParser(description="Do compartment analysis for group of cool files")
parser.add_argument("--cools", nargs="+", required=True, help="cool/mcool files")
parser.add_argument("--tags", nargs="+", help="Names for the cool files")
parser.add_argument("--bin", type=int, default=20000)
parser.add_argument("--region", default="/mnt/data0/wenliang/database/mm10/mm10.sizes.bed")
parser.add_argument("--outdir", default=".", help="Directory to put the results in")
args = parser.parse_args()

chromsizes = bioframe.fetch_chromsizes('mm10')
# chromsizes = pd.read_csv(args.size, sep="\t", index_col=0, names=["start", "end"])
chromosomes = list(chromsizes.index)
chrom_regions = pd.concat([pd.Series([0] * len(chromosomes), index=chromsizes.index, name='start'), chromsizes], axis=1).reset_index()

if args.tags:
	tags = args.tags
else:
	tags = []
	for coolfile in args.cools:
		filename = os.path.basename(coolfile)
		tags.append(re.sub("\.\S+", "", filename))

conditions = tags
binsize = args.bin

bin_coolfile = [f'{x}::/resolutions/{binsize}' if x.endswith("mcool") else x for x in args.cools]

cooler_paths = dict(zip(tags, bin_coolfile))
long_names = dict(zip(tags, tags))
pal = sns.color_palette('colorblind')
colors = {tags[i] : pal[i] for i in range(len(tags))}

clrs = {
	cond: cooler.Cooler(cooler_paths[cond]) for cond in conditions
}


f, axs = plt.subplots(figsize=(8,13), nrows=2, gridspec_kw={'height_ratios':[6,2]}, sharex=True)

for cond in conditions:
	clr = clrs[cond]
	# chromsizes = chromsizes.set_index('chrom')
	cvd = cooltools.expected.diagsum(
		clr = clrs[cond],
		supports = chromsizes,
		transforms={'balanced': lambda p: p['count']*p['weight1']*p['weight2']}
	)

	cvd_agg = (
		cvd.groupby('diag').agg({'n_valid':'sum', 'count.sum':'sum', 'balanced.sum':'sum',}).reset_index()
	)

	lb_cvd, lb_slopes, lb_distbins = cooltools.expected.logbin_expected(cvd)
	lb_cvd_agg, lb_slopes_agg = cooltools.expected.combine_binned_expected(lb_cvd, binned_exp_slope=lb_slopes)
	lb_cvd_agg['s_bp'] = lb_cvd_agg['diag.avg'] * clr.binsize 
	lb_slopes_agg['s_bp'] = lb_slopes_agg['diag.avg'] * clr.binsize

	cvd_agg['s_bp'] = (cvd_agg['diag'] * clr.binsize)
	cvd_agg['count.avg'] = (cvd_agg['count.sum'] / cvd_agg['n_valid'])
	cvd_agg['balanced.avg'] = (cvd_agg['balanced.sum'] / cvd_agg['n_valid'])

	ax = axs[0]
	ax.loglog(lb_cvd_agg['s_bp'], lb_cvd_agg['balanced.avg'], 'o-', markersize=5, color=colors[cond])

	ax.set(ylabel='IC contact frequency', xlim=(1e3, 1e8))
	ax.set_aspect(1.0)
	ax.grid(lw=0.5)

	ax = axs[1]
	ax.semilogx(lb_slopes_agg['s_bp'], lb_slopes_agg['slope'], alpha=0.5, color=colors[cond])
	ax.set(xlabel='separation, bp', ylabel='slope')
	ax.grid(lw=0.5)

plt.savefig(f"{args.outdir}/contacts_vs_distance.pdf")
