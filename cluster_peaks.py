import os
import re
import argparse
import pybedtools


def consec(lst, dist, length):
	it = iter(lst)
	prev = next(it)
	tmp = [prev]
	for ele in it:
		if prev[0] != ele[0] or int(ele[1]) - int(prev[2]) > dist or int(ele[2]) - int(tmp[0][1]) > length:
			yield tmp
			tmp = [ele]
		else:
			tmp.append(ele)
		prev = ele
	yield tmp


def cluster_beds(bedfile, dist, length, outdir, prefix):
	lst = pybedtools.BedTool(bedfile)
	lst = lst.sort()

	clusters = consec(lst, dist, length)
	if prefix:
		filename = prefix
	else:
		filename = os.path.basename(bedfile)
		filename = re.sub("\.\S+", "", filename)
	out_fh = open(f"{outdir}/{filename}.cluster.bed", 'w')
	for cluster in clusters:
		size = len(cluster)
		chrom = cluster[0][0]
		start = cluster[0][1]
		end = cluster[-1][2]
		print(f'{chrom}\t{start}\t{end}\t{size}', file=out_fh)


def main():
	parser = argparse.ArgumentParser(description="Cluster bed regions")
	parser.add_argument("bedfile", help="Bedfile that have the regions")
	parser.add_argument("-d", type=int, default=20000, help="Distance between regions")
	parser.add_argument("-l", type=int, default=200000, help="lenght of the cluster")
	parser.add_argument("-o", default=".", help="output dir")
	parser.add_argument("-p", help="file name")
	args = parser.parse_args()

	cluster_beds(args.bedfile, args.d, args.l, args.o, args.p)


if __name__ == '__main__':
	main()
