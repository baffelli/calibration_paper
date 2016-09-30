import matplotlib.pyplot as plt
import pyrat.core.corefun as cf
import pyrat.fileutils.gpri_files as gpf
import pyrat.visualization.visfun as vf


def plot_figure_1(inputs, outputs, threads, config, params, wildcards):
    plt.style.use(inputs['style'])
    # This is a full page figure, so we create a figure twice as wide as usual
    # create figure
    fig_w, fig_h = plt.rcParams['figure.figsize']
    f, axarr = plt.subplots(2, 3, figsize=(fig_w * 3, fig_h * 2))
    for idx_chan, chan_str in enumerate(['HH', 'VV']):
        for idx_proc, proc_str in enumerate(['', '_desq', '_corr']):
            current_ax = axarr[idx_chan, idx_proc]
            file_name = inputs[chan_str + proc_str]
            current_slc = gpf.gammaDataset(file_name + '.par', file_name)  # load slc
            ptarg_zoom, r_plot, az_plot, mx_pos, res_dict, r_vec, az_vec = cf.ptarg(current_slc, params['ridx'],
                                                                                    params['azidx'],
                                                                                    azwin=128, rwin=20, osf=32)
            mph, rgb, norm = vf.dismph(ptarg_zoom, k=0.3, sf=0.2, coherence=False)  # create rgb image
            mappable = current_ax.imshow(mph, cmap=rgb, extent=[az_vec.min(), az_vec.max(), r_vec.min(), r_vec.max()],
                                         aspect=1 / 10.0)
            # add resolution analysis
            resolution_text = r"""Range resolution: {r_res:.3f} m
            Azimuth resolution: {az_res:.3f}$^\circ$""".format(r_res=res_dict['range_resolution'][0],
                                             az_res=res_dict['azimuth_resolution'][
                                                 0])
            # print(resolution_text)
            bbox_props = dict(boxstyle="square", fc="white", ec="w", lw=2)
            t = current_ax.text(0.1, 0.1, resolution_text, size=8, bbox=bbox_props, horizontalalignment='left',
                                transform=current_ax.transAxes)  # set label
            if idx_chan == 1 and idx_proc == 0:
                current_ax.set_xlabel(r"azimuth samples $\theta$ [deg]")
                current_ax.set_ylabel(r'range samples [m]')

    ticks = [0, 0.25, 0.5, 0.75, 1]
    labels = ['-180', '-90', '0', '90', '180']
    # Make space for colorbar
    # f.subplots_adjust(right=0.8)
    # cbar_ax = divider.append_axes("left", size="20%", pad=0.25)
    cax = f.add_axes([0.25, 0.02, 0.5, 0.015])
    cbar = f.colorbar(mappable, label=r'Phase [deg]', orientation='horizontal', ticks=ticks, cax=cax)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(labels)
    cbar.update_ticks()
    f.savefig(outputs[0])


plot_figure_1(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)
