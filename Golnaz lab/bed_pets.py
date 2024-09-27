import os
import re
import numpy as np
import pandas as pd
import cooler
import argparse


def count_gene_pets(coolfile, gene_df, extend=5000):
	# group_gene_df = gene_df.groupby('chrom')
	cool = cooler.Cooler(coolfile)
	chromosomes = ['chr' + str(x) for x in range(1, 19)]
	chromosomes = chromosomes + ["chrX", "chrY"]

	gene_df = gene_df[gene_df['chrom'].isin(chromosomes)]
	group_gene_df = gene_df.groupby('chrom')

	m, n = gene_df.shape
	gene_pet_list = [0] * m
	for chrom, chrom_gene_df in group_gene_df:
		chrom_matrix = cool.matrix(balance=False, as_pixels=True, join=True).fetch(chrom)

		for i, row in chrom_gene_df.iterrows():
			up_gene_matrix = chrom_matrix.loc[(chrom_matrix['start1'] >= row.start - extend) & \
				(chrom_matrix['end1'] <= row.end + extend)]
			down_gene_matrix = chrom_matrix.loc[(chrom_matrix['start2'] >= row.start - extend) & \
				(chrom_matrix['end2'] <= row.end + extend)]

			gene_matrix = pd.concat([up_gene_matrix, down_gene_matrix])
			gene_matrix = gene_matrix.drop_duplicates()

			gene_pets = gene_matrix['count'].sum()
			index = gene_df.index[gene_df['geneID'] == row.geneID].tolist()
			gene_pet_list[index[0]] = gene_pets

	sample_name = re.sub("\.\S+", "", os.path.basename(coolfile))
	gene_df[sample_name] = gene_pet_list
	return gene_df


def cal_logFC(row):
	try:
		log2FC = np.log2(row.TCF1_3T3 / row.Norm_3T3)
	except:
		log2FC = "NA"

	return log2FC


def main():
	parser = argparse.ArgumentParser(description="Calculate PETs in gene body")
	parser.add_argument("cool", nargs="+", help="cool files")
	parser.add_argument("--gene", help="gene annotation file")
	parser.add_argument("--extend", default=5000, type=int)
	parser.add_argument("--outdir", default=".")
	args = parser.parse_args()

	gene_df = pd.read_csv(args.gene, sep="\t", names=['chrom', 'start', 'end', 'geneID', 'geneName'])

	for coolfile in args.cool:
		gene_df = count_gene_pets(coolfile, gene_df, extend=args.extend)

	# print(gene_df.head())
	gene_df['logFC'] = gene_df.apply(lambda row: cal_logFC(row), axis=1)
	gene_df.to_csv(os.path.join(args.outdir, "gene_pets.csv"), sep="\t", index=False)


if __name__ == '__main__':
	main()
