import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
import pyrat.core.corefun as cf
import pyrat.diff.intfun as _ifun
import pyrat.fileutils.gpri_files as gpf
import pyrat.visualization.visfun as vf
from matplotlib import gridspec

from matplotlib.patches import Ellipse

# Percentile
perc = 20
# Coherence threshold
coh_thresh = 0.6

win = [5, 3]  # multilooking window


# Return the decimated azimuth position
def az_idx(ds, idx):
    return idx / ds.GPRI_decimation_factor


# Confidence for the error ellipse
conf = 5.991

# SLice for image
sl_image = (slice(500 // win[0], 1900 // win[0]), slice(None, None))


def plot_figure_9(inputs, outputs, threads, config, params, wildcards):
    HHVV = gpf.gammaDataset(inputs.C_cal_par, inputs.HHVV, dtype=gpf.type_mapping['FCOMPLEX'])
    HH = gpf.gammaDataset(inputs.C_cal_par, inputs.HH, dtype=gpf.type_mapping['FCOMPLEX'])
    VV = gpf.gammaDataset(inputs.C_cal_par, inputs.VV, dtype=gpf.type_mapping['FCOMPLEX'])
    # height
    theta = gpf.load_binary(inputs.theta, VV.shape[0], dtype=gpf.type_mapping['FLOAT'])
    # EStimate coherence
    HHVV = _ifun.estimate_coherence(HHVV, HH, VV, win, discard=True)
    mli = cf.smooth(VV, win, discard=True)
    theta = cf.smooth(theta, win, discard=True)
    # Copolar span
    copol_span = mli
    # Take the brightest 10%
    bright_percentile = np.percentile(copol_span, perc)
    # Find the indices of that percentile
    perc_r, perc_az = np.nonzero((copol_span > bright_percentile) * (np.abs(HHVV) > coh_thresh) * (theta != 0))
    # Create RGB
    # Property of the mph image
    mph_dict = {'k': 0.08, 'sf': 0.9, 'coherence': True, 'peak': False, 'mli': mli[sl_image], 'coherence_threshold': 0.4,
                'coherence_slope': 6}
    HHVV_im = HHVV[sl_image]
    az_vec = HHVV_im.az_vec
    r_vec = HHVV_im.r_vec
    mph, rgb, norm = vf.dismph(HHVV_im, **mph_dict)  # create rgb image
    pal, ext = vf.dismph_palette(HHVV_im, **mph_dict)
    # Set style and create figure
    plt.style.use(inputs['style'])
    fig_w, fig_h = plt.rcParams['figure.figsize']
    f = plt.figure(figsize=(fig_w*2, fig_h*2))
    # Create grid of plots
    gs = gridspec.GridSpec(*(3, 2), height_ratios=[1, 1, 0.2])
    gs.update(hspace=0.3, wspace=0.1)
    im_ax = f.add_subplot(gs[0:2, ::])
    aspect = fig_h / fig_w
    slc_ext = [az_vec[0], az_vec[-1], r_vec[-1], r_vec[1]]
    # Show the phase
    im_ax.imshow(mph, extent=slc_ext, aspect=vf.fixed_aspect(slc_ext, aspect), origin='upper')
    # Show the indices of bright targets
    # im_ax.scatter(az_vec[perc_az], r_vec[perc_r], facecolors='none', edgecolors='cyan')
    im_ax.yaxis.set_label_text(r'Range [m]')
    im_ax.xaxis.set_label_text(r'Azimuth angle [$^\circ$]')
    im_ax.yaxis.set_major_locator(tick.MultipleLocator(500))
    im_ax.xaxis.set_major_locator(tick.MultipleLocator(20))
    im_ax.set_title('HH-VV Phase')
    # # Plot reflectors
    # pos_list = {'Simmleremoos 2': (-20, 15), 'Simmleremoos 1': (-15, -15)}  # position to avoid overlapping
    box = dict(boxstyle="round", fc="w", lw=0.2)
    # for ref in params['ref']:
    #     dec_pos = (int(ref['ridx']) / win[0], HHVV.azidx_dec(int(ref['azidx'])))
    #     grid_pos = (az_vec[dec_pos[1]], r_vec[dec_pos[0]])
    #     im_ax.plot(*grid_pos, marker='o', markeredgecolor='#feb24c', mfc='none', mew=1, ms=10)
    #     annotations = plt.annotate(ref['name'], xy=grid_pos, color='black', size=7,
    #                                xytext=pos_list.get(ref['name'], (0, -15)), textcoords='offset points', bbox=box,
    #                                horizontalalignment='center')
    # plot palette
    pal_ax = f.add_subplot(gs[-1, 0])
    pal_ax.imshow(pal, aspect=1 / 30.0, extent=[0, 1, -np.pi, np.pi, ])
    pal_ax.set_ylabel(r'Phase')
    pal_ax.set_xlabel(r'Intensity')
    pal_ax.grid(b=False)
    pal_ax.set_yticks([-np.pi, 0, np.pi])
    pal_ax.set_yticklabels([r"$-\pi$", r"$0$", r"$\pi$"])
    ext_list = [ext[0], (ext[0] + ext[1]) / 2, ext[1]]
    pal_ax.set_xticks([])
    # pal_ax.set_xticklabels([r"{{val:.2f}}".format(val=val) for val in ext_list])
    # Plot coherence
    coh_ax = f.add_subplot(gs[-1,1])
    c = np.linspace(0, 1)
    c_scale = vf.scale_coherence(c, threshold=mph_dict['coherence_threshold'], slope=mph_dict['coherence_slope'])
    coh_ax.plot(c, c_scale)
    coh_ax.set_yticks([0, 0.5, 1])
    coh_ax.xaxis.set_label_text(r'Coherence')
    coh_ax.yaxis.set_label_text(r'Saturation')
    coh_ax.set_aspect(1)


    # Plot correlation of height with copol phase
    hist_fig, hist_ax = plt.subplots(figsize=(fig_w, fig_h))
    HHVV_bright = np.angle(HHVV[perc_r, perc_az])
    theta_bright = np.sin(theta[perc_r, perc_az].real)
    coherence_bright = np.abs(HHVV[perc_r, perc_az])
    # Compute covariance
    # #Normalie
    theta_var = np.std(theta_bright)
    HHVV_var = np.std(HHVV_bright)
    theta_bright_norm = (theta_bright - np.mean(theta_bright)) / theta_var
    HHVV_bright_norm = (HHVV_bright - np.mean(HHVV_bright)) / HHVV_var
    x_cov = np.cov(theta_bright_norm, y=HHVV_bright_norm)
    # # Ellispe parameters
    # if wildcards.n == '10':
    #     l, w = np.linalg.eigh(x_cov)
    #     pos = [0.4, 0.6]
    #     hist_ax.annotate("correlation: {:1.1f}".format(x_cov[1,0]), xy=pos, bbox=box, xytext=pos,
    #                      textcoords='axes fraction', xycoords='axes fraction', fontsize=plt.rcParams['axes.labelsize'])
    hist_ax.hist2d(theta_bright, HHVV_bright, bins=200)
    hist_ax.set_xlabel(r' $\sin\left(\theta_l\right)$')
    hist_ax.set_ylabel(r'HH-VV phase [rad]')
    # hist_ax.set_aspect(1)
    # hist_ax.set_aspect('square')
    # ticks = [-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi]
    # hist_ax.set_yticks(ticks)
    # hist_ax.set_yticklabels([r"$-\pi$", r"$-\frac{\pi}{2}$", r"$0$", r"$\frac{\pi}{2}$", r"$\pi$"])
    f.savefig(outputs['fig_a'])
    hist_fig.savefig(outputs['fig_b'])


plot_figure_9(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)
