import os
import sys
import re
import pandas as pd
import numpy as np
import networkx as nx
import pybedtools
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import argparse


def find_community(bedpefile, outdir, prefix):
	exits_nodes = dict()
	G = nx.Graph()
	for line in open(bedpefile, 'r'):
		row = line.rstrip().split()
		row[1] = int(row[1])
		row[2] = int(row[2])
		row[4] = int(row[4])
		row[5] = int(row[5])

		nodeA = f'{row[0]}\t{row[1]}\t{row[2]}'
		exits_nodes_str = "\n".join(exits_nodes.keys())
		intersect = pybedtools.BedTool(nodeA, from_string=True).intersect(pybedtools.BedTool(exits_nodes_str, from_string=True), wo=True)
		if intersect:
			for hit in intersect:
				nodeA = f'{hit[3]}\t{hit[4]}\t{hit[5]}'

		nodeB = f'{row[3]}\t{row[4]}\t{row[5]}'
		intersect = pybedtools.BedTool(nodeB, from_string=True).intersect(pybedtools.BedTool(exits_nodes_str, from_string=True), wo=True)
		if intersect:
			for hit in intersect:
				nodeB = f'{hit[3]}\t{hit[4]}\t{hit[5]}'

		G.add_edge(nodeA, nodeB, weight=1)
		exits_nodes[nodeA] = 1
		exits_nodes[nodeB] = 1

	communities = list(reversed(sorted(list(nx.connected_components(G)), key=len)))
	sizes = [len(x) for x in communities]
	rank = range(1, len(sizes) + 1)
	node_weight = dict(G.degree())

	outfile = os.path.join(outdir, f"{prefix}.Anchor_info.bed")
	outfh = open(outfile, 'w')
	i = 0
	for community in communities:
		i += 1
		size = len(community)
		for node in community:
			node_w = node_weight[node]
			print(f'{node}\t{node_w}\t{size}\tC{i}', file=outfh)

	ax = plt.subplot(111)
	ax.scatter(rank, sizes, s=0.3)
	# ax.set_xlim([0, 100])
	# ax.set_ylim([0, 50])
	ax.set_xlabel("Rank")
	ax.set_ylabel("Size of community")
	plt.savefig(f"{outdir}/{prefix}.size_rank.pdf")


def main():
	parser = argparse.ArgumentParser(description="Detect the 3D community from Hi-C loops")
	parser.add_argument("bedpe", help="The bedpe file that have the loops")
	parser.add_argument("--outdir", default=".")
	parser.add_argument("--prefix")
	args = parser.parse_args()

	if not args.prefix:
		prefix = re.sub(".bedpe", "", os.path.basename(args.bedpe))
	else:
		prefix = args.prefix

	find_community(args.bedpe, args.outdir, prefix)


if __name__ == '__main__':
	main()
