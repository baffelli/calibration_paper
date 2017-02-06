import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import pyrat.core.corefun as cf
import pyrat.fileutils.gpri_files as gpf

# Window for plotting in degress
sw_az = 2
sw_az_plot = 0.8
sw_r = 2


def plot_figure_2(inputs, outputs, threads, config, params, wildcards):
    plt.style.use(inputs['style'])
    slc_HH = gpf.gammaDataset(inputs['VV'][0] + '.par', inputs['VV'][0])
    slc_VV = gpf.gammaDataset(inputs['VV'][1] + '.par', inputs['VV'][1])
    fig_w, fig_h = plt.rcParams['figure.figsize']
    f, phase_ax = plt.subplots(1, sharex=True, figsize=(fig_w, fig_h))
    sorted_by_range = sorted(params['reflectors'], key=lambda tup: tup['ridx'])
    # sorted_by_range = [ref for ref in sorted_by_range if ref[] == "t"]
    cm = mpl.cm.get_cmap('inferno', len(sorted_by_range))  # colormap for sorting
    ridx_vec = [el['ridx'] for el in sorted_by_range]
    ranges = slc_VV.r_vec[ridx_vec]
    bd = np.hstack((ranges, ranges[
        -1] * 2))  # stupid trick, the last boundary is twice the last range to that we can set the tik to be in the middle of the cell
    nm = mpl.colors.BoundaryNorm(bd, cm.N)
    mappable = mpl.cm.ScalarMappable(cmap=cm, norm=nm)
    mappable.set_array(bd)  # set the array of range distances to the mappable
    line = {}
    for slc, chan_name, ls in zip((slc_VV, slc_HH), ['VV', 'HH'], ['-', '--']):
        line[chan_name] = []
        for ref_vec in sorted_by_range:
            sw_az_id = sw_az / slc.GPRI_az_angle_step  # number of samples in degrees
            sw_az_id_plot = sw_az_plot / slc.GPRI_az_angle_step  # number of samples in degrees
            # Extract maximum around the supposed position
            ridx_max, azidx_max = cf.maximum_around(np.abs(slc), [ref_vec['ridx'], ref_vec['azidx']], [sw_r, sw_az_id])
            # Slice the slc
            slc_sl = (ridx_max, slice(azidx_max - sw_az_id_plot / 2, azidx_max + sw_az_id_plot / 2))
            reflector_slice = slc[slc_sl]
            max_slc_VV = slc[ridx_max, azidx_max]
            # Azimuth angle vector for plot
            az_vec = slc.GPRI_az_angle_step * np.arange(-len(reflector_slice) / 2
                                                        , len(reflector_slice) / 2)
            refl_ph = np.angle(reflector_slice) - np.angle(max_slc_VV)
            refl_ph = np.unwrap(refl_ph)
            refl_amp = (np.abs(reflector_slice)) ** 2
            r_sl = slc.r_vec[ref_vec['ridx']]
            line[chan_name].append(phase_ax.plot(az_vec, np.rad2deg(refl_ph), ls=ls, color=mappable.to_rgba(r_sl))[0])
            # Only format VV
            lab = r"{name}, {chan}".format(name=ref_vec['name'], chan=chan_name)
            # amp_ax.plot(az_vec, refl_amp / refl_amp[reflector_slice.shape[0] / 2], color=mappable.to_rgba(r_sl),
            #             label=lab, ls=ls)
    # Add dummies for legend

    # Plot line for beamwidth
    phase_ax.set_ylim(-30, 30)
    lc = "#cc8e8e"
    lw = 1
    phase_ax.axvline(0.2, color=lc, ls='-', lw=lw)
    phase_ax.axvline(-0.2, color=lc, ls='-', lw=lw)
    phase_ax.yaxis.set_label_text(r'Phase [deg]')
    phase_ax.axvline(0.2, color=lc, ls='-', lw=lw)
    phase_ax.axvline(-0.2, color=lc, ls='-', lw=lw)
    # ph.yaxis.set_label_text(r'Intensity')
    phase_ax.xaxis.set_label_text(r'Azimuth angle from maximum [$\circ$]')
    # Format axes
    # map(vf.format_axes, f.get_axes())
    # leg = phase_ax.legend([line[0][0], line[0][1]] + line[:][0], ['HH', 'VV'] + [ref['name'] for ref in sorted_by_range],
    #                       ncol=len(ranges) // 2)
    lines_leg = [(l_h, l_v) for l_h,l_v in zip(line['HH'], line['VV'])]
    # labels =[ref['name'] for ref in sorted_by_range]
    # leg_line = plt.legend(lines_leg, labels, ncol=len(ranges) // 2, loc='lower center')
    # leg_pol = plt.legend( [line['HH'][0] ,line['VV'][0]] ,['HH','VV'])
    # phase_ax.add_artist(leg_line)
    # phase_ax.add_artist(leg_pol)
    phase_ax.legend([line['HH'][0], line['VV'][0]] + lines_leg,
                    ['HH', 'VV'] + [ref['name'] for ref in sorted_by_range], ncol=2, fancybox=True, frameon=True, shadow=False,
                          framealpha=None, numpoints=2)
    # leg = phase_ax.legend(loc='lower center', ncol=len(ranges) // 2,
    # leg.get_title().set_fontsize(plt.rcParams['legend.fontsize'])
    # cbar_ax = f.add_axes([0.95, 0.15, 0.03, 0.7])
    # cbar = f.colorbar(mappable=mappable, cax=cbar_ax, boundaries=bd, norm=nm, orientation='horizontal')#colorbar encoding distance
    # tick position is at the midpoint of ranges
    # tick_pos = (bd[1:] + bd[:-1])/2
    # cbar.set_ticks(tick_pos)
    # cbar.set_ticklabels(["{:.2f}".format(r) for r in ranges])
    # cbar.set_clim(ranges[0], ranges[-1])
    # cbar.set_label('Distance from radar [m]')
    # plt.show()
    f.subplots_adjust(top=1)
    f.savefig(outputs[0])
    plt.close(f)


plot_figure_2(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)
