import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
import pyrat.core.corefun as cf
import pyrat.diff.intfun as _ifun
import pyrat.fileutils.gpri_files as gpf
import pyrat.visualization.visfun as vf
from matplotlib import gridspec

import scipy.stats as stats

# Percentile
perc = 20
# Coherence threshold
coh_thresh = 0.6

win = [3, 2]  # multilooking window

# Return the decimated azimuth position
def az_idx(ds, idx):
    return idx / ds.GPRI_decimation_factor


def plot_figure_9(inputs, outputs, threads, config, params, wildcards):
    print(inputs.HHVV_phase)
    HHVV = gpf.gammaDataset(inputs.C_cal_par, inputs.HHVV_phase[0], dtype=gpf.type_mapping['FCOMPLEX'])
    HH = gpf.gammaDataset(inputs.C_cal_par, inputs.HHVV_phase[1], dtype=gpf.type_mapping['FCOMPLEX'])
    VV = gpf.gammaDataset(inputs.C_cal_par, inputs.HHVV_phase[2], dtype=gpf.type_mapping['FCOMPLEX'])
    #height
    theta = np.pi/2 + gpf.load_binary(inputs.theta, VV.shape[0], dtype=gpf.type_mapping['FLOAT'])
    # u = gpf.load_binary(inputs.u, VV.shape[0], dtype=gpf.type_mapping['FLOAT'])
    # theta = theta - u
    #Topographic phase
    topo_phase = gpf.gammaDataset(inputs.topo_phase_par, inputs.topo_phase, dtype=gpf.type_mapping['FCOMPLEX'])
    topo_phase = cf.smooth(topo_phase, win, discard=True)
    # EStimate coherence
    HHVV = _ifun.estimate_coherence(HHVV, HH, VV, win, discard=True)
    mli = cf.smooth(VV, win, discard=True)
    theta = cf.smooth(theta, win, discard=True)
    print(HHVV.shape)
    print(topo_phase.shape)
    az_vec = HHVV.az_vec
    r_vec = HHVV.r_vec
    # Copolar span
    copol_span = mli
    # Take the brightest 10%
    bright_percentile = np.percentile(copol_span, perc)
    # Find the indices of that percentile
    perc_r, perc_az = np.nonzero((copol_span > bright_percentile) * (np.abs(HHVV) > coh_thresh))
    # perc_subs = np.unravel_index(perc_indices, copol_span.shape)

    # Create RGB
    mph_dict = {'k': 0.08, 'sf': 0.9, 'coherence': True, 'peak': False, 'mli': mli, 'coherence_threshold': 0.6,
                'coherence_slope': 12}
    mph, rgb, norm = vf.dismph(HHVV, **mph_dict)  # create rgb image
    pal, ext = vf.dismph_palette(HHVV, **mph_dict)
    plt.style.use(inputs['style'])
    fig_w, fig_h = plt.rcParams['figure.figsize']
    f = plt.figure(figsize=(2 * fig_w, 2 * fig_h))
    gs = gridspec.GridSpec(*(2, 4), height_ratios=[1, 0.2])
    gs.update(hspace=0.3, wspace=0.5)
    im_ax = f.add_subplot(gs[0, 2::])
    aspect = fig_h / fig_w
    slc_ext = [az_vec[0], az_vec[-1], r_vec[-1], r_vec[1]]
    # Azimuth range grid
    aa, rr = np.meshgrid(az_vec, r_vec, indexing='xy')
    # Show the phase
    im_ax.imshow(mph, extent=slc_ext, aspect=vf.fixed_aspect(slc_ext, aspect), origin='upper')
    # Show the indices of bright targets
    # im_ax.scatter(az_vec[perc_az], r_vec[perc_r], facecolors='none', edgecolors='cyan')
    im_ax.yaxis.set_label_text(r'Range [m]')
    im_ax.xaxis.set_label_text(r'Azimuth angle [$^\circ$]')
    im_ax.yaxis.set_major_locator(tick.MultipleLocator(500))
    im_ax.xaxis.set_major_locator(tick.MultipleLocator(20))
    im_ax.set_title('HH-VV Phase')
    # Plot reflectors
    pos_list = {'Simmleremoos 2': (-20, 15), 'Simmleremoos 1': (-15, -15)}  # position to avoid overlapping
    box = dict(boxstyle="round", fc="w", lw=0.2)
    for ref in params['ref']:
        dec_pos = (int(ref['ridx']) / win[0], HHVV.azidx_dec(int(ref['azidx'])))
        grid_pos = (az_vec[dec_pos[1]], r_vec[dec_pos[0]])
        im_ax.plot(*grid_pos, marker='o', markeredgecolor='#feb24c', mfc='none', mew=1, ms=10)
        annotations = plt.annotate(ref['name'], xy=grid_pos, color='black', size=7,
                                   xytext=pos_list.get(ref['name'], (0, -15)), textcoords='offset points', bbox=box,
                                   horizontalalignment='center')
    # plot palette
    pal_ax = f.add_subplot(gs[-1, 2])
    print(ext)
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
    coh_ax = f.add_subplot(gs[-1, 3])
    c = np.linspace(0, 1)
    c_scale = vf.scale_coherence(c, threshold=mph_dict['coherence_threshold'], slope=mph_dict['coherence_slope'])
    coh_ax.plot(c, c_scale)
    coh_ax.xaxis.set_label_text(r'Coherence')
    coh_ax.yaxis.set_label_text(r'Saturation')
    coh_ax.set_aspect(1)
    #Plot correlation of height with copol phase
    hist_ax = f.add_subplot(gs[0, 1])
    HHVV_bright = np.angle(HHVV[perc_r, perc_az])
    topo_bright = np.angle(topo_phase[perc_r, perc_az])
    theta_bright = theta[perc_r, perc_az].real
    coherence_bright = np.abs(HHVV[perc_r, perc_az])
    hist_ax.hist2d(theta_bright,HHVV_bright, bins=100)
    hist_ax.set_xlabel(r'Elevation angle')
    hist_ax.set_ylabel(r'HH-VV Phase')
    #Fit
    slope, intercept, r,p, err = stats.linregress(theta_bright.flatten(),y=HHVV_bright.flatten())
    print(r)
    hist_ax.plot(theta_bright, slope*theta_bright+intercept,np.pi)
    plt.show()

    # # Plot Histogram of Copol Phase
    # hist_ax = f.add_subplot(gs[0, 1])
    # HHVV_bright = np.angle(HHVV[perc_r, perc_az])
    # HHVV_mean = np.mean(HHVV_bright)
    # #PLot histogram
    # hist_ax.hist(HHVV_bright, 100, range=(-np.pi, np.pi), stacked=True, normed=True)
    # #Plot average
    # hist_ax.axvline(x=HHVV_mean,color='g')
    # hist_ax.set_ylim([0, 1.5])
    # ticks = [-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi]
    # hist_ax.set_xticks(ticks)
    # hist_ax.set_xticklabels([r"$-\pi$", r"$-\frac{\pi}{2}$", r"$0$", r"$\frac{\pi}{2}$", r"$\pi$"])
    # hist_ext = [-np.pi, np.pi, 0, 1.5]
    # hist_ax.set_aspect(vf.fixed_aspect(hist_ext, aspect))
    # # hist_ax.set_yticks([])
    # # hist_ax.set_yticklabels([])
    # hist_ax.set_xlabel(r'$arg\left(\gamma_{HHVV}\right)$')
    # hist_ax.set_ylabel(r'Normalized frequency')
    # title_str = r'HH-VV phase on brightest {perc:1.1f}\%  pixels with $\vert\gamma_{{HHVV}}\vert $\ge$ {coh_thresh}'.format(
    #     perc=100 - perc, coh_thresh=coh_thresh)
    # hist_ax.set_title(title_str)
    f.savefig(outputs[0])


plot_figure_9(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)
