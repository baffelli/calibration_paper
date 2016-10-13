import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pyrat.core.corefun as cf
import pyrat.fileutils.gpri_files as gpf
import pyrat.visualization.visfun as vf

# Window for plotting in degress
ws = 1.2
sw = [2, 3]  # search window for maxium


def plot_figure_2(inputs, outputs, threads, config, params, wildcards):
    plt.style.use(inputs['style'])
    slc = gpf.gammaDataset(inputs['VV'] + '.par', inputs['VV'])
    fig_w, fig_h = plt.rcParams['figure.figsize']
    f, (phase_ax, amp_ax) = plt.subplots(2, sharex=True, figsize=(fig_w * 2, fig_h))
    sorted_by_range = sorted(params['reflectors'], key=lambda tup: tup['ridx'])
    # sorted_by_range = [ref for ref in sorted_by_range if ref[] == "t"]
    cm = mpl.cm.get_cmap('viridis', len(sorted_by_range))  # colormap for sorting
    ridx_vec = [el['ridx'] for el in sorted_by_range]
    ranges = slc.r_vec[ridx_vec]
    nm = mpl.colors.Normalize(vmin=ranges[0], vmax=ranges[-1])
    mappable = mpl.cm.ScalarMappable(cmap=cm, norm=nm)
    mappable.set_array(ranges)  # set the array of range distances to the mappable
    for ref_vec in sorted_by_range:

        ws_id = ws / slc.GPRI_az_angle_step[0]  # number of samples in degrees
        # Extract maximum around the supposed position
        ridx_max, azidx_max = cf.maximum_around(np.abs(slc), [ref_vec['ridx'], ref_vec['azidx']], [sw[0], ws_id])
        # Slice the slc
        slc_sl = (ridx_max, slice(azidx_max - ws_id / 2, azidx_max + ws_id / 2))
        reflector_slice = slc[slc_sl]
        max_slc = slc[ridx_max, azidx_max]
        # Azimuth angle vector for plot
        az_vec = slc.GPRI_az_angle_step[0] * np.arange(-len(reflector_slice) / 2
                                                       , len(reflector_slice) / 2)
        refl_ph = np.angle(reflector_slice) - np.angle(max_slc)
        refl_ph = np.unwrap(refl_ph)
        refl_amp = (np.abs(reflector_slice)) ** 2
        r_sl = slc.r_vec[ref_vec['ridx']]
        line, = phase_ax.plot(az_vec, np.rad2deg(refl_ph), color=mappable.to_rgba(r_sl))
        amp_ax.plot(az_vec, refl_amp / refl_amp[reflector_slice.shape[0] / 2], color=mappable.to_rgba(r_sl),
                    label=r"{r_sl:.0f}".format(r_sl=r_sl))
        # Plot line for beamwidth
        phase_ax.set_ylim(-30, 30)
        lc = "#cc8e8e"
        lw = 0.5
        phase_ax.axvline(0.2, color=lc, ls='-', lw=lw)
        phase_ax.axvline(-0.2, color=lc, ls='-', lw=lw)
        phase_ax.yaxis.set_label_text(r'Phase [deg]')
        amp_ax.axvline(0.2, color=lc, ls='-', lw=lw)
        amp_ax.axvline(-0.2, color=lc, ls='-', lw=lw)
        amp_ax.yaxis.set_label_text(r'Intensity')
        amp_ax.xaxis.set_label_text(r'azimuth angle from maximum [deg]')
        # Format axes
        map(vf.format_axes, f.get_axes())
    leg = amp_ax.legend(loc='lower center', ncol=len(ranges) // 2, fancybox=True, frameon=True, shadow=False, framealpha=None,
                  title=r"Distance from radar [m]", handlelength=0.3)
    leg.get_title().set_fontsize(plt.rcParams['legend.fontsize'])
    f.subplots_adjust(top=1)
    f.savefig(outputs[0])
    plt.close(f)


plot_figure_2(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)
