import os
import re
import argparse
from subprocess import check_output
import pandas as pd


def domain_binding(domain, bamfile, outdir):
	binding_strength = list()

	bedfile = os.path.join(outdir, "binding_strength.bed")
	domain_fh = open(domain, 'r')
	domain_fh.readline()
	for line in domain_fh.readlines():
		row = line.rstrip().split()
		row[1] = int(row[1])
		row[2] = int(row[2])
		length = row[2] - row[1]
		cmd = 'bamCoverage -r {chrom}:{start}:{end} -bs 10000 --normalizeUsing RPKM -b {bam} -of bedgraph -o {bed}'.format(chrom=row[0], \
			start=row[1], end=row[2], bam=bamfile, bed=bedfile)
		cmd_out = check_output(cmd, shell=True)
		
		bed_df = pd.read_csv(bedfile, sep="\t", names=['chrom', 'start', 'end', 'rpkm'])
		RPKM = bed_df['rpkm'].mean()

		binding_strength.append(RPKM)
	domain_fh.close()
	return binding_strength

def domain_peaks(domain, peakfile):
	peak_density = list()
	cmd = 'bedtools intersect -a {} -b {} -wao'.format(domain, peakfile)
	cmd_out = check_output(cmd, shell=True, universal_newlines=True)
	last_pos_key = ""
	last_peak_num = 0
	seen_peak = dict()
	# peak_pattern = re.compile(r"[b-z]$")
	for line in cmd_out.rstrip().split("\n"):
		row = line.rstrip().split()
		pos_key = "{}\t{}\t{}".format(row[0], row[1], row[2])
		# row[8] = re.sub("[a-z]$", "", row[8])
		row[-1] = int(row[-1])
		# if peak_pattern.match(row[8]):
		#	continue
		
		if last_pos_key != pos_key and last_pos_key != "":
			chrom, start, end = last_pos_key.split("\t")
			start = int(start)
			end = int(end)
			bin_size = end - start
			density = 1000000*float(last_peak_num)/float(bin_size)
			peak_density.append(density)
			last_peak_num = 0
		
		if row[-1] > 0:
			last_peak_num += 1
		last_pos_key = pos_key

	chrom, start, end = last_pos_key.split("\t")
	start = int(start)
	end = int(end)
	bin_size = end - start
	density = float(last_peak_num)/float(bin_size)
	peak_density.append(density)

	return peak_density

def main():
	parser = argparse.ArgumentParser(description="Calculate the binding strength of protein within TADs")
	parser.add_argument("domain", help="File that have the domain information")
	# parser.add_argument("bam", help="Bam file of ChIP seq")
	parser.add_argument("peak", help="Peak file of ChIP seq")
	parser.add_argument("--outdir", default=".")
	parser.add_argument("--key", default="Peak")
	args = parser.parse_args()

	domain_df = pd.read_csv(args.domain, sep="\t")
	# binding_strength = domain_binding(args.domain, args.bam, args.outdir)
	# domain_df['RPKM'] = binding_strength

	peak_density = domain_peaks(args.domain, args.peak)
	domain_df[args.key] = peak_density

	domain_df.to_csv(os.path.join(args.outdir, "Domain_score_binding.tab"), sep="\t", index=False)

if __name__ == '__main__':
	main()
