import os
import re
import pandas as pd
import numpy as np
from collections import Counter, defaultdict
import argparse


def consec(lst, dist, length):
	it = iter(lst)
	prev = next(it)
	tmp = [prev]
	for ele in it:
		if prev + dist <= ele or ele - tmp[0] > length:
			yield tmp
			tmp = [ele]
		else:
			tmp.append(ele)
		prev = ele
	yield tmp


def count_peaks(peakfile, chrom, start, end):
	pass


def integrate_bound(insulations, tags):
	is_boundary = defaultdict(list)
	coordinates = defaultdict(list)
	strength = defaultdict(list)
	for i, insulation in enumerate(insulations):
		cell = tags[i]
		ins_df = pd.read_csv(insulation, sep="\t")
		boundary_df = ins_df.loc[ins_df['boundary_strength_500000']>=0.1]
		for x in boundary_df.index:
			is_boundary[x].append(cell)
			coordinates[x] = [ins_df.iloc[x]['chrom'], ins_df.iloc[x]['start'], ins_df.iloc[x]['end']]
			strength[x].append(ins_df.iloc[x]['boundary_strength_500000'])

	merged_i = list(consec(sorted(list(is_boundary.keys())), dist=5, length=10))

	# outfh = open(outfile, 'w')
	stat = defaultdict(dict)
	for block in merged_i:
		if len(block) > 1:
			block_key = f'{block[0]}-{block[-1]}'
			coords = f'{coordinates[block[0]][0]}\t{coordinates[block[0]][1]}\t{coordinates[block[-1]][2]}'
			cells = [j for x in block for j in is_boundary[x]]
			cell_str = ",".join(cells)
			block_str = ",".join(map(str, block))
			boundary_strength = np.mean([j for x in block for j in strength[x]])
			# print(f"M{block_key}\t{coords}\t{cell_str}\t{boundary_strength}", file=outfh)

			stat[block_key]['chrom'] = coordinates[block[0]][0]
			stat[block_key]['start'] = coordinates[block[0]][1]
			stat[block_key]['end'] = coordinates[block[-1]][2]
			stat[block_key]['strength'] = boundary_strength
			for cell in cells:
				stat[block_key][cell] = 1
		else:
			block_key = f'{block[0]}'
			cell_str = ",".join(is_boundary[block[0]])
			coords = f'{coordinates[block[0]][0]}\t{coordinates[block[0]][1]}\t{coordinates[block[0]][2]}'
			boundary_strength = np.mean([j for j in strength[block[0]]])
			# print(f"{block_key}\t{coords}\t{cell_str}\t{boundary_strength}", file=outfh)

			stat[block_key]['chrom'] = coordinates[block[0]][0]
			stat[block_key]['start'] = coordinates[block[0]][1]
			stat[block_key]['end'] = coordinates[block[0]][2]
			stat[block_key]['strength'] = boundary_strength
			for cell in is_boundary[block[0]]:
				stat[block_key][cell] = 1
	stat_df = pd.DataFrame.from_dict(stat, orient='index')
	stat_df = stat_df.fillna(0)
	stat_df.index.name = 'boundaryID'
	stat_df.to_csv("boundary.stat.csv", index=False)

def main():
	parser = argparse.ArgumentParser(description="parse insulation files and integrate boundaries")
	parser.add_argument("-ins", nargs="+", required=True)
	parser.add_argument("-tags", nargs="+")
	# parser.add_argument("-out", default='integrated_bound.tab')
	args = parser.parse_args()

	tags = []
	if args.tags:
		tags = args.tags
	else:
		for file in args.ins:
			filename = os.path.basename(file)
			tag = re.sub("\.\S+", "", filename)
			tags.append(tag)

	integrate_bound(args.ins, tags)


if __name__ == '__main__':
	main()
