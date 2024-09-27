import os
import re
import cooler
import pandas as pd
import numpy as np
import straw
import argparse


def Compute_Score_hic(tadfile, hicfile, res):
	domain_scores = list()
	domain_pets = list()
	tad_fh = open(tadfile, 'r')
	for line in tad_fh:
		row = line.rstrip().split()
		row[1] = int(row[1])
		row[2] = int(row[2])
		chrom_i = re.sub("chr", "", row[0])
		pos_string = "{}:{}:{}".format(chrom_i, row[1], row[2])

		dump_hic = straw.straw("NONE", hicfile, pos_string, chrom_i, "BP", res)
		matrix_df = pd.DataFrame(np.array(dump_hic).T)
		matrix_df.columns = ["start", "end", "contact"]
		matrix_df = matrix_df.astype(int)

		intraTAD_df = matrix_df.loc[(matrix_df['start'] >= row[1]) & (matrix_df['end'] <= row[2])]

		all_contacts = matrix_df['contact'].sum()
		intraTAD_contacts = intraTAD_df['contact'].sum()

		try:
			score = intraTAD_contacts / all_contacts
		except:
			score = "NA"
		domain_scores.append(score)
		domain_pets.append(intraTAD_contacts)

	return domain_scores, domain_pets


def Compute_Score_cool(tadfile, coolfile):
	cool = cooler.Cooler(coolfile)
	tad_df = pd.read_csv(tadfile, sep="\t", names=['Chrom', 'Start', 'End'])
	group_tad_df = tad_df.groupby("Chrom")
	m, n = tad_df.shape
	domain_scores = [0] * m
	domain_pets = [0] * m
	for chrom, chrom_tad_df in group_tad_df:
		chrom_matrix = cool.matrix(balance=False, as_pixels=True, join=True).fetch(chrom)

		for i, row in chrom_tad_df.iterrows():
			up_matrix_df = chrom_matrix.loc[(chrom_matrix['start1'] >= row.Start) & (chrom_matrix['start1'] <= row.End)]
			down_matrix_df = chrom_matrix.loc[(chrom_matrix['start2'] >= row.Start) & (chrom_matrix['start2'] <= row.End)]
			matrix_df = pd.concat([up_matrix_df, down_matrix_df])
			matrix_df = matrix_df.drop_duplicates()

			intraTAD_df = matrix_df.loc[(matrix_df['start1'] >= row.Start) & (matrix_df['start2'] <= row.End)]

			all_contacts = matrix_df['count'].sum()
			intraTAD_contacts = intraTAD_df['count'].sum()

			try:
				score = float(intraTAD_contacts) / float(all_contacts)
			except:
				score = "NA"
			index = tad_df.index[(tad_df['Chrom'] == chrom) & (tad_df['Start'] == row.Start)].tolist()
			domain_scores[index[0]] = score
			domain_pets[index[0]] = intraTAD_contacts

	return domain_scores, domain_pets


def main():
	parser = argparse.ArgumentParser(description="Compute the domain score for a common set of TADs")
	parser.add_argument("tads", help="common TADs in bed format")
	parser.add_argument("--hic", nargs="+", help="HiC files for different conditions")
	parser.add_argument("--cool", nargs="+", help="Cool files for different conditions")
	parser.add_argument("--res", type=int, default=10000, help="resolution")
	parser.add_argument("--outdir", default=".", help="output directory")
	parser.add_argument("--juicer", default="/mnt/data0/wenliang/software/juicer/scripts/juicer_tools.jar", help="path for juice_tools")
	args = parser.parse_args()

	domain_scores_df = pd.read_csv(args.tads, sep="\t", names=["Chrom", "Start", "End"])
	domain_pets_df = pd.read_csv(args.tads, sep="\t", names=['Chrom', 'Start', 'End'])

	if args.hic:
		for hicfile in args.hic:
			tag = re.sub("\.\S+", "", os.path.basename(hicfile))
			domain_scores, domain_pets = Compute_Score_hic(args.tads, hicfile, res=10000)
			domain_scores_df[tag] = domain_scores
			domain_pets_df[tag] = domain_pets

	if args.cool:
		for coolfile in args.cool:
			filename = os.path.basename(coolfile)
			tag = re.sub("\.\S+", "", filename)

			if filename.endswith(".cool"):
				domain_scores, domain_pets = Compute_Score_cool(args.tads, coolfile)
			elif filename.endswith(".mcool"):
				coolfile = "{}::/resolutions/{}/".format(coolfile, args.res)
				domain_scores, domain_pets = Compute_Score_cool(args.tads, coolfile)
			domain_scores_df[tag] = domain_scores
			domain_pets_df[tag] = domain_pets

	# domain_scores_df = pd.DataFrame.from_dict(domain_scores)
	domain_scores_df.to_csv(os.path.join(args.outdir, "Domain_scores.tab"), sep="\t", index=False)
	domain_pets_df.to_csv(os.path.join(args.outdir, "Domain_pets.tab"), sep="\t", index=False)


if __name__ == '__main__':
	main()
