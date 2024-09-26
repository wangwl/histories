import numpy as np
import os
import re
import argparse


def get_intra(matrix):
	n_rows, n_cols = matrix.shape

	up_left = matrix[0: n_rows // 2, 0: n_rows // 2]
	lower_right = matrix[(n_rows // 2) + 1: n_rows, (n_rows // 2) + 1: n_rows]

	up_mean = np.nanmean(up_left)
	lower_mean = np.nanmean(lower_right)

	avg_mean = (up_mean + lower_mean) / 2.0
	return avg_mean


def get_inter(matrix):
	n_rows, n_cols = matrix.shape
	# upright_matrix = matrix[(2 * n_rows// 3) +1 : n_rows, 0 : n_rows// 3]
	upright_matrix = matrix[(n_rows // 2) + 1: n_rows, 0: n_rows // 2]
	# print(upright_matrix.shape)
	mean = np.nanmean(upright_matrix)
	return mean


def main():
	parser = argparse.ArgumentParser(description="Calculate the enrichment of upright square")
	parser.add_argument("--matrix", required=True, nargs="+")
	args = parser.parse_args()

	for matrixfile in args.matrix:
		matrix = np.loadtxt(matrixfile)
		inter = round(get_inter(matrix), 3)
		intra = round(get_intra(matrix), 3)
		row = os.path.basename(matrixfile).split(".")
		cell = row[0]
		cond = ".".join(row[1:-1])
		print("{}\t{}\t{}\t{}".format(cell, cond, inter, intra))


if __name__ == '__main__':
	main()
