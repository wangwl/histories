import os
import argparse


def get_TADs(insulation):
	insulation_fh = open(insulation, 'r')
	insulation_fh.readline()

	select_chroms = ['chr' + str(x) for x in range(1,20)]
	select_chroms.append("chrX")
	last_chrom = ""
	TAD_st = None
	TAD_ed = None
	for line in insulation_fh:
		row = line.rstrip().split()

		if len(row) == 0:
			continue

		if row[0] != last_chrom and TAD_st != None:
			print("{}\t{}\t{}".format(last_chrom, TAD_st, TAD_ed))
			TAD_st = None

		if row[0] not in select_chroms:
			continue

		if TAD_st == None and row[3] =='False':
			TAD_st = int(row[1])

		if float(row[6]) > 0.1:
			TAD_ed = int(row[1])
			print("{}\t{}\t{}".format(row[0], TAD_st, TAD_ed))
			TAD_st = int(row[2])

		last_chrom = row[0]
		TAD_ed = int(row[2])

def main():
	parser = argparse.ArgumentParser(description="Generate TADs from the insulation file")
	parser.add_argument("insulation")
	args = parser.parse_args()

	get_TADs(args.insulation)


if __name__ == '__main__':
	main()
