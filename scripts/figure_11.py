import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
import pyrat.core.corefun as cf
import pyrat.diff.intfun as _ifun
import pyrat.fileutils.gpri_files as gpf
import pyrat.visualization.visfun as vf
from matplotlib import gridspec
import matplotlib.ticker as tick
import string
import scipy.signal as _sig


# Return the decimated azimuth position
def az_idx(ds, idx):
    return idx / ds.azimuth_looks

def center_width_from_slice(sl):
    center = (sl.stop + sl.start)//2
    width = (sl.stop - sl.start)
    return center, width

k = 0.2
sf = 0.7
mph_dict = {'k': k, 'sf': sf, 'coherence': False, 'peak': False}
az_sl = slice(1392,3600)
r_sl = slice(550,1850)



def plot_figure_11(inputs, outputs, threads, config, params, wildcards):
    # Create figure
    plt.style.use(inputs['style'])
    fig_w, fig_h = plt.rcParams['figure.figsize']
    aspect = fig_h / fig_w
    f = plt.figure(figsize=(2 * fig_w, 2 * fig_h))
    # Grid of plots for figures
    gs = gridspec.GridSpec(*(2, 3), height_ratios=[1,0.2])
    gs.update(hspace=0.0, wspace=0.25)
    #data to plot
    slc_names = ['slc', 'slc_desq', 'slc_corr']
    for plot_idx, slc_name in enumerate(slc_names):
        #Load slc
        slc = gpf.gammaDataset(inputs[slc_name + '_par'], inputs[slc_name])
        az_sl_dec = slice(az_sl.start//slc.GPRI_decimation_factor, az_sl.stop//slc.GPRI_decimation_factor)
        slc = slc[r_sl, az_sl_dec]
        slc_ext = [slc.az_vec[0], slc.az_vec[-1], slc.r_vec[-1], slc.r_vec[0]]
        slc_aspect = slc_ext[1] - slc_ext[0]
        #Plot raw channel
        # raw_ax = f.add_subplot(gs[0,plot_idx])
        # raw_ax.imshow(vf.exp_im(np.abs(_sig.hilbert(raw[raw_sl],axis=0)),k,sf), extent=raw_ext, aspect=1e-5, interpolation='none', cmap='gray')
        #Plot slc
        mph, rgb, norm = vf.dismph(slc, **mph_dict)  # create rgb image
        pal, ext = vf.dismph_palette(slc, **mph_dict)
        slc_ax = f.add_subplot(gs[::, plot_idx])
        print(vf.fixed_aspect(slc_ext, aspect))
        slc_ax.imshow(mph, extent=slc_ext, aspect= vf.fixed_aspect(slc_ext, aspect),  origin='upper', interpolation='none', cmap='gray')
        slc_ax.xaxis.set_major_locator(tick.MultipleLocator(10))
        tit = string.ascii_lowercase[plot_idx]
        slc_ax.title.set_text(r"({label_name})".format(label_name=tit))
        if plot_idx == 0:
            slc_ax.set_xlabel(r"Azimuth [$^\circ$]")
            slc_ax.set_ylabel(r'Range [m]')
    pal_ax = f.add_subplot(gs[-1,1])
    pal_ax.imshow(pal, extent=ext, aspect=2)
    pal_ax.set_ylabel(r'Phase')
    pal_ax.set_xlabel(r'Intensity')
    pal_ax.grid(b=False)
    pal_ax.set_yticks([-np.pi, 0, np.pi])
    pal_ax.set_yticklabels([r"$-\pi$", r"$0$", r"$\pi$"])
    ext_list = [ext[0], (ext[0] + ext[1]) / 2, ext[1]]
    pal_ax.set_xticks(ext_list)
    f.savefig(outputs[0])
    # plt.show()
    # print(inputs.HHVV_phase)
    # HHVV = gpf.gammaDataset(inputs.C_cal_par, inputs.HHVV_phase[0], dtype=gpf.type_mapping['FCOMPLEX'])
    # HH = gpf.gammaDataset(inputs.C_cal_par, inputs.HHVV_phase[1], dtype=gpf.type_mapping['FCOMPLEX'])
    # VV = gpf.gammaDataset(inputs.C_cal_par, inputs.HHVV_phase[2], dtype=gpf.type_mapping['FCOMPLEX'])
    # win = [3, 2]  # multilooking window
    # HHVV = _ifun.estimate_coherence(HHVV, HH, VV, win, discard=True)
    # mli = cf.smooth(VV, win, discard=True)
    # az_vec = HHVV.az_vec
    # r_vec = HHVV.r_vec
    # # Create RGB
    # mph_dict = {'k': 0.08, 'sf': 0.9, 'coherence': True, 'peak': False, 'mli': mli, 'coherence_threshold': 0.6,
    #             'coherence_slope': 12}
    # mph, rgb, norm = vf.dismph(HHVV, **mph_dict)  # create rgb image
    # pal, ext = vf.dismph_palette(HHVV, **mph_dict)
    # print(ext)
    #
    # plt.style.use(inputs['style'])
    # fig_w, fig_h = plt.rcParams['figure.figsize']
    # f = plt.figure(figsize=(2 * fig_w, 2 * fig_h))
    # gs = gridspec.GridSpec(*(2, 4), height_ratios=[1, 0.2])
    # gs.update(hspace=0.3, wspace=0.0)
    # im_ax = f.add_subplot(gs[0, ::])
    # aspect = fig_h / fig_w
    # im_ax.imshow(mph, extent=[az_vec[0], az_vec[-1], r_vec[-1], r_vec[1]], aspect=1 / 20 * aspect, origin='upper')
    # im_ax.yaxis.set_label_text(r'range [m]')
    # im_ax.xaxis.set_label_text(r'azimuth [$^\circ$]')
    # im_ax.yaxis.set_major_locator(tick.MultipleLocator(500))
    # im_ax.xaxis.set_major_locator(tick.MultipleLocator(20))
    # # Plot reflectors
    # pos_list = {'Simmleremoos 2': (-10, 15), 'Simmleremoos 1': (-10, -15)}  # position to avoid overlapping
    # box = dict(boxstyle="round", fc="w", lw=0.2)
    # for ref in params['ref']:
    #     dec_pos = (int(ref['ridx']) / win[0], HHVV.azidx_dec(int(ref['azidx'])))
    #     grid_pos = (az_vec[dec_pos[1]], r_vec[dec_pos[0]])
    #     im_ax.plot(*grid_pos, marker='o', markeredgecolor='#43a2ca', mfc='none', mew=2, ms=10)
    #     annotations = plt.annotate(ref['name'], xy=grid_pos, color='black', size=7,
    #                                xytext=pos_list.get(ref['name'], (0 ,-15)), textcoords='offset points', bbox=box,
    #                                horizontalalignment='center')
    # # plot palette
    # pal_ax = f.add_subplot(gs[-1, 1])
    # pal_ax.imshow(pal, extent=ext, aspect=500)
    # pal_ax.set_ylabel(r'Phase')
    # pal_ax.set_xlabel(r'Intensity')
    # pal_ax.grid(b=False)
    # pal_ax.set_yticks([-np.pi, 0, np.pi])
    # pal_ax.set_yticklabels([r"$-\pi$", r"$0$", r"$\pi$"])
    # ext_list = [ext[0], (ext[0] + ext[1]) / 2, ext[1]]
    # pal_ax.set_xticks([])
    # # pal_ax.set_xticklabels([r"{{val:.2f}}".format(val=val) for val in ext_list])
    # # Plot coherence
    # coh_ax = f.add_subplot(gs[-1, 2])
    # c = np.linspace(0, 1)
    # c_scale = vf.scale_coherence(c, threshold=mph_dict['coherence_threshold'], slope=mph_dict['coherence_slope'])
    # coh_ax.plot(c, c_scale)
    # coh_ax.xaxis.set_label_text(r'Coherence')
    # coh_ax.yaxis.set_label_text(r'Saturation')
    # coh_ax.set_aspect(1)
    # f.savefig(outputs[0])


plot_figure_11(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)
