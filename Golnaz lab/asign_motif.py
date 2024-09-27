import bioframe
import os
import argparse
import subprocess
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict


JASPAR_URL = "http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2020/mm10/"

def asign_motif(peakfile, motif):
	peaks = bioframe.read_table(peakfile, schema='bed5')
	jaspar_motif_file = f'{motif}.tsv.gz'
	motifs = bioframe.read_table(JASPAR_URL + jaspar_motif_file, schema='jaspar', skiprows=1)

	peaks_motifs = bioframe.overlap(peaks, motifs)
	peaks_motifs['pval_2'] = peaks_motifs['pval_2'].fillna(-1) 
	idxmax_peaks_motifs = peaks_motifs.groupby(["chrom_1", "start_1","end_1"])["pval_2"].idxmax().values
	df_peaks_maxmotif = peaks_motifs.loc[idxmax_peaks_motifs]
	df_peaks_maxmotif['pval_2'].replace(-1,pd.NA,inplace=True)

	df_peaks_maxmotif = df_peaks_maxmotif[['chrom_1', 'start_1', 'end_1', 'start_2', 'end_2', 'pval_2', 'strand_2']]
	df_peaks_maxmotif.columns = ['chrom', 'start', 'end', 'mstart', 'mend', 'pval', 'strand']

	return df_peaks_maxmotif


def overlap_peaks(peak_motif_df, peak2_motif2_df):
	overlap_peaks = bioframe.overlap(peak_motif_df, peak2_motif2_df)
	
	distance = list()
	strands = list()
	strand_distance = defaultdict(list)
	for i, row in overlap_peaks.iterrows():
		try:
			motif_coords = sorted([row.mstart_1, row.mend_1, row.mstart_2, row.mend_2])
			dist = motif_coords[3] - motif_coords[0]
			strand = row.strand_1 + row.strand_2
		except:
			dist = 0
			strand = 'no'
		distance.append(dist)
		strands.append(strand)
		strand_distance[strand].append(dist)

	motif_dist = pd.DataFrame(list(zip(distance, strands)), columns=['Distance', 'Strands'])
	# motif_dist = pd.DataFrame.from_dict(strand_distance, orient='columns')
	motif_dist.to_csv("motif_dist.csv", index=False)
	sns.kdeplot(motif_dist, x="Distance", hue='Strands', multiple="stack")
	plt.savefig("Motif.dist.pdf")


def main():
	parser = argparse.ArgumentParser(description="Find the motif of each peak")
	parser.add_argument("peakfile", help="Bed file that has all the peaks")
	parser.add_argument("motif", help="The JASPAR motif ID")
	parser.add_argument("--peak2", help="Another peak file")
	parser.add_argument("--motif2", help="Motif ID of the second TF")
	args = parser.parse_args()

	peak_motif_df = asign_motif(args.peakfile, args.motif)

	if args.peak2 and args.motif2:
		peak2_motif2_df = asign_motif(args.peak2, args.motif2)
		overlap_peaks(peak_motif_df, peak2_motif2_df)


if __name__ == '__main__':
	main()
		