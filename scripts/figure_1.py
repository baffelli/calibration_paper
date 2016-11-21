import string

import matplotlib.pyplot as plt
import numpy as np
import pyrat.core.corefun as cf
import pyrat.fileutils.gpri_files as gpf
import pyrat.visualization.visfun as vf
from matplotlib import gridspec

# searc window in degrees
sw = 0.8
# window for interpolation
rwin = 20
azwin = 128


def plot_figure_1(inputs, outputs, threads, config, params, wildcards):
    plt.style.use(inputs['style'])
    # This is a full page figure, so we create a figure twice as wide as usual
    # create figure
    fig_w, fig_h = plt.rcParams['figure.figsize']
    rs = 2
    shp = (2 * rs + 1, 3)#shape of gridspec
    f = plt.figure(figsize=(fig_w * 3, fig_h * 2))
    gs = gridspec.GridSpec(*shp)
    gs.update(wspace=0.0, hspace=0.7)
    # f, ax_arr = plt.subplots(nrows=3, ncols=3, figsize=(fig_w * 3, fig_h * 2), gridspec_kw={'height_ratios':[1,1,0.5]})
    for idx_chan, chan_str in enumerate(['HH', 'VV']):
        for idx_proc, proc_str in enumerate(['', '_desq', '_corr']):
            # current label
            start_idx = idx_chan + (idx_chan%rs)
            current_ax = plt.subplot(gs[start_idx:start_idx+rs, idx_proc])
            current_idx = np.ravel_multi_index((idx_chan , idx_proc), shp)
            label_name = string.ascii_lowercase[current_idx]
            # current_ax = axarr[idx_chan, idx_proc]
            current_ax.title.set_text(r"({label_name})".format(label_name=label_name))
            file_name = inputs[chan_str + proc_str]
            current_slc = gpf.gammaDataset(file_name + '.par', file_name)  # load slc
            sw_idx = sw / current_slc.GPRI_az_angle_step
            ptarg_zoom, r_plot, az_plot, mx_pos, res_dict, r_vec, az_vec = cf.ptarg(current_slc, params['ridx'],
                                                                                    params['azidx'],
                                                                                    azwin=azwin, rwin=rwin, osf=32,
                                                                                    sw=(4, sw_idx))
            # remove phase at maximum
            ptarg_zoom *= np.exp(1j * np.angle(ptarg_zoom[mx_pos].conj()))
            mph_dict = {'k':0.35, 'sf':0.3, 'coherence':False, 'peak':True}
            mph, rgb, norm = vf.dismph(ptarg_zoom, **mph_dict)  # create rgb image
            pal, ext = vf.dismph_palette(ptarg_zoom, **mph_dict)
            aspect = fig_h / fig_w
            new_aspect = rwin / azwin * aspect  # in this way the subplot has the same aspect as the figure
            mappable = current_ax.imshow(mph, cmap=rgb, extent=[az_vec.min(), az_vec.max(), r_vec.min(), r_vec.max()],
                                         aspect=new_aspect)
            # add resolution analysis
            resolution_text = r"""Range resolution: {r_res:.3f} m
            Azimuth resolution: {az_res:.3f}$^\circ$""".format(r_res=res_dict['range_resolution'][0],
                                                               az_res=res_dict['azimuth_resolution'][
                                                                   0])
            bbox_props = dict(boxstyle="square", fc="white", ec="w", lw=2)
            t = current_ax.text(0.1, 0.1, resolution_text, size=8, bbox=bbox_props, horizontalalignment='left',
                                transform=current_ax.transAxes)  # set label
            if idx_chan == 1 and idx_proc == 0:
                current_ax.set_xlabel(r"Azimuth [deg]")
                current_ax.set_ylabel(r'Range [m]')

    # Make space for colorbar
    # f.subplots_adjust(hspace=0.3, wspace=0.2)
    # cax = f.add_axes([0.5, 0, 0.5, 0.07], anchor=(0.5,0.5))
    cax = plt.subplot(gs[-1,:])
    cax.imshow(pal, aspect=1/30.0, extent=[ 0, 1, -np.pi, np.pi,])
    cax.set_ylabel(r'Phase')
    cax.set_xlabel(r'Intensity relative to peak')
    cax.grid(b=False)
    cax.set_yticks([-np.pi, 0, np.pi])
    cax.set_yticklabels([r"$-\pi$", r"$0$",  r"$\pi$"])
    cax.set_xticks([0, 0.5, 1])
    cax.set_xticklabels([r"$0$", r"$0.5$",  r"$1$"])

    map(vf.format_axes, f.get_axes())
    f.savefig(outputs[0])


plot_figure_1(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)
