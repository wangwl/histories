import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
import fanc
import fanc.plotting
from scipy import ndimage as ndi
import matplotlib.patches as patches
from scipy.ndimage import zoom
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import re
import argparse
from subprocess import check_output


def clipped_zoom(img, zoom_factor, **kwargs):
	h, w = img.shape[:2]
	zoom_tuple = (zoom_factor,) * 2 + (1,) * (img.ndim - 2)
	if zoom_factor < 1:
		zh = int(np.round(h * zoom_factor))
		zw = int(np.round(w * zoom_factor))
		top = (h - zh) // 2
		left = (w - zw) // 2
		out = np.zeros_like(img)
		out[top: top + zh, left: left + zw] = zoom(img, zoom_tuple, **kwargs)
	elif zoom_factor > 1:
		zh = int(np.round(h / zoom_factor))
		zw = int(np.round(w / zoom_factor))
		top = (h - zh) // 2
		left = (w - zw) // 2
		out = zoom(img[top: top + zh, left: left + zw], zoom_tuple, **kwargs)
		trim_top = ((out.shape[0] - h) // 2)
		trim_left = ((out.shape[1] - w) // 2)
		out = out[trim_top: trim_top + h, trim_left: trim_left + w]
	else:
		out = img
	return out


def highlight_features(dataframe, region, color, a, axes):
	try:
		features = dataframe.loc[region].values.tolist()
		if type(features[0]) == int:
			_, x_min, x_max, y_min, y_max = features
			rect = patches.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min, linewidth=1.2, edgecolor=color, facecolor='none')
			axes[a].add_patch(rect)
		else:
			for f in features:
				_, x_min, x_max, y_min, y_max = f
				rect = patches.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min, linewidth=1.2, edgecolor=color, facecolor='none')
				axes[a].add_patch(rect)
	except KeyError:
		next


def visualize(chess, pairs, hicA, hicB, outdir, features, sn_thr, zsim_thr):
	similarities = pd.read_csv(chess, sep='\t', index_col=0)
	regions = pd.read_csv(pairs, sep='\t', header=None)
	sub_sim = similarities[(similarities["SN"] >= sn_thr) & (similarities['z_ssim'] <= zsim_thr)]

	patient_hic = fanc.load(hicA)
	control_hic = fanc.load(hicB)
	patient_tag = re.sub("\.\S+", "", os.path.basename(hicA))
	control_tag = re.sub("\.\S+", "", os.path.basename(hicB))

	gain_feature = os.path.join(features, "gained_features.tsv")
	loss_feature = os.path.join(features, "lost_features.tsv")
	if not os.path.exists(gain_feature):
		regions2compare = regions.loc[sub_sim.index]
		regions2compare.to_csv('filtered_{}'.format(pairs), '\t', index=False, header=False)

		filtered_file = 'filtered_{}'.format(pairs)
		extract_cmd = f'chess extract {filtered_file} {hicA} {hicB} {features}'
		check_output(extract_cmd, shell=True)
	gained = pd.read_csv(gain_feature, delimiter=',', usecols=[0, 1, 2, 3, 4, 5], header=None, index_col=[0], engine='python')
	lost = pd.read_csv(loss_feature, delimiter=',', usecols=[0, 1, 2, 3, 4, 5], header=None, index_col=[0], engine='python')

	cluster_ids = list()
	cluster_ids.append(sub_sim.index[0])

	plot_dir = f'{outdir}/plots'
	if not os.path.exists(plot_dir):
		os.mkdir(plot_dir)
	for region_id in sub_sim.index[1:]:
		if region_id - 1 == cluster_ids[-1]:
			cluster_ids.append(region_id)
			continue
		else:
			chrom, window_start = regions.loc[cluster_ids[0]][0:2]
			window_end = regions.loc[cluster_ids[-1]][2]
			region_string = "{}:{}-{}".format(chrom, window_start, window_end)
			patient_region_sub = patient_hic[region_string, region_string].data
			control_region_sub = control_hic[region_string, region_string].data
			min_v = min([np.min(np.extract(patient_region_sub > 0, patient_region_sub)), np.min(np.extract(control_region_sub > 0, control_region_sub))])
			patient_region_sub += min_v
			control_region_sub += min_v
			l2fcm = np.log2(patient_region_sub / control_region_sub)
			zml2 = clipped_zoom(l2fcm, 0.7)
			rot_l2 = ndi.rotate(zml2, 45, reshape=False)
			fig, axes = plt.subplots(1, 3, figsize=(9, 3))
			axes[0].set_title(patient_tag)
			axes[1].set_title(control_tag)
			axes[2].set_title(f'$log_2$({patient_tag} / {control_tag})')
			zm1 = clipped_zoom(control_region_sub, 0.7)
			rot_control = ndi.rotate(zm1, 45, reshape=False)
			zm2 = clipped_zoom(patient_region_sub, 0.7)
			rot_patient = ndi.rotate(zm2, 45, reshape=False)
			middle = int(np.shape(rot_control)[1] / 2.)
			m1 = axes[0].imshow(rot_patient[:middle, :], vmin=0, vmax=0.03, cmap='germany')
			m2 = axes[1].imshow(rot_control[:middle, :], vmin=0, vmax=0.03, cmap='germany')
			# for reg in cluster_ids:
			# 	highlight_features(gained, reg, 'crimson', 0, axes)
			# 	highlight_features(lost, reg, 'royalblue', 1, axes)
			m3 = axes[2].imshow(rot_l2[:middle, :], cmap='seismic', vmax=5, vmin=-5)
			for m, ax in zip([m1, m2, m3], axes):
				ax.axis('off')
				divider = make_axes_locatable(ax)
				cax = divider.append_axes('bottom', size='5%', pad=0.05)
				fig.colorbar(m, cax=cax, orientation='horizontal')
			plt.tight_layout()
			plt.savefig("{}/{}_{}_{}.pdf".format(plot_dir, chrom, window_start, window_end))
			cluster_ids = [region_id]

	chrom, window_start = regions.loc[cluster_ids[0]][0:2]
	window_end = regions.loc[cluster_ids[-1]][2]
	region_string = "{}:{}-{}".format(chrom, window_start, window_end)
	patient_region_sub = patient_hic[region_string, region_string].data
	control_region_sub = control_hic[region_string, region_string].data
	min_v = min([np.min(np.extract(patient_region_sub > 0, patient_region_sub)), np.min(np.extract(control_region_sub > 0, control_region_sub))])
	patient_region_sub += min_v
	control_region_sub += min_v
	l2fcm = np.log2(patient_region_sub / control_region_sub)
	zml2 = clipped_zoom(l2fcm, 0.7)
	rot_l2 = ndi.rotate(zml2, 45, reshape=False)
	fig, axes = plt.subplots(1, 3, figsize=(9, 3))
	axes[0].set_title(patient_tag)
	axes[1].set_title(control_tag)
	axes[2].set_title(f'$log_2$({patient_tag} / {control_tag})')
	zm1 = clipped_zoom(control_region_sub, 0.7)
	rot_control = ndi.rotate(zm1, 45, reshape=False)
	zm2 = clipped_zoom(patient_region_sub, 0.7)
	rot_patient = ndi.rotate(zm2, 45, reshape=False)
	middle = int(np.shape(rot_control)[1] / 2.)
	m1 = axes[0].imshow(rot_patient[:middle, :], vmin=0, vmax=0.03, cmap='germany')
	m2 = axes[1].imshow(rot_control[:middle, :], vmin=0, vmax=0.03, cmap='germany')
	for reg in cluster_ids:
		highlight_features(gained, reg, 'crimson', 0, axes)
		highlight_features(lost, reg, 'royalblue', 1, axes)
	m3 = axes[2].imshow(rot_l2[:middle, :], cmap='seismic', vmax=5, vmin=-5)
	for m, ax in zip([m1, m2, m3], axes):
		ax.axis('off')
		divider = make_axes_locatable(ax)
		cax = divider.append_axes('bottom', size='5%', pad=0.05)
		fig.colorbar(m, cax=cax, orientation='horizontal')
	plt.tight_layout()
	plt.savefig("{}/{}_{}_{}.pdf".format(plot_dir, chrom, window_start, window_end))


def main():
	parser = argparse.ArgumentParser(description="visualize the diff regions from chess")
	parser.add_argument("chess", help="Chess sim results")
	parser.add_argument("pairs", help="Pair file from chess pairs")
	parser.add_argument("hicA", help="The first hic file given to chess sim")
	parser.add_argument("hicB", help="The second hic file given to chess sim")
	parser.add_argument("--outdir", default=".", help="output direcotry")
	parser.add_argument("--features", help="the features direcotry")
	parser.add_argument("--sn", type=float, default=0.5)
	parser.add_argument("--zsim", type=float, default=-1)
	args = parser.parse_args()

	if not args.features:
		features = os.path.join(args.outdir, 'features')
	else:
		features = args.features

	visualize(args.chess, args.pairs, args.hicA, args.hicB, args.outdir, features, args.sn, args.zsim)


if __name__ == '__main__':
	main()
