import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pyrat.core.corefun as cf
import pyrat.fileutils.gpri_files as gpf
import pyrat.visualization.visfun as vf

# Window for plotting in degress
sw_az = 2
sw_az_plot = 0.8
sw_r = 2

def plot_figure_2(inputs, outputs, threads, config, params, wildcards):
    plt.style.use(inputs['style'])
    slc_HH = gpf.gammaDataset(inputs['VV'][0] + '.par', inputs['VV'][0])
    slc_VV = gpf.gammaDataset(inputs['VV'][1] + '.par', inputs['VV'][1])
    fig_w, fig_h = plt.rcParams['figure.figsize']
    f, (phase_ax, amp_ax) = plt.subplots(2, sharex=True, figsize=(fig_w * 2, fig_h * 2))
    sorted_by_range = sorted(params['reflectors'], key=lambda tup: tup['ridx'])
    # sorted_by_range = [ref for ref in sorted_by_range if ref[] == "t"]
    cm = mpl.cm.get_cmap('inferno', len(sorted_by_range))  # colormap for sorting
    ridx_vec = [el['ridx'] for el in sorted_by_range]
    ranges = slc_VV.r_vec[ridx_vec]
    bd = np.hstack((ranges, ranges[-1] *2))#stupid trick, the last boundary is twice the last range to that we can set the tik to be in the middle of the cell
    nm = mpl.colors.BoundaryNorm(bd, cm.N)
    mappable = mpl.cm.ScalarMappable(cmap=cm, norm=nm)
    mappable.set_array(bd)  # set the array of range distances to the mappable
    for ref_vec in sorted_by_range:


        for slc, chan_name, ls in zip((slc_VV, slc_HH), ['VV', 'HH'], ['-', '--']):
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
            line, = phase_ax.plot(az_vec, np.rad2deg(refl_ph), color=mappable.to_rgba(r_sl), ls=ls)
            #Only format VV
            lab = r"{name}, {chan}".format(name=ref_vec['name'], chan=chan_name)
            amp_ax.plot(az_vec, refl_amp / refl_amp[reflector_slice.shape[0] / 2], color=mappable.to_rgba(r_sl),
                        label=lab, ls=ls)
        # Plot line for beamwidth
        phase_ax.set_ylim(-30, 30)
        lc = "#cc8e8e"
        lw = 1
        phase_ax.axvline(0.2, color=lc, ls='-', lw=lw)
        phase_ax.axvline(-0.2, color=lc, ls='-', lw=lw)
        phase_ax.yaxis.set_label_text(r'Phase [deg]')
        amp_ax.axvline(0.2, color=lc, ls='-', lw=lw)
        amp_ax.axvline(-0.2, color=lc, ls='-', lw=lw)
        amp_ax.yaxis.set_label_text(r'Intensity')
        amp_ax.xaxis.set_label_text(r'azimuth angle from maximum [deg]')
        # Format axes
        map(vf.format_axes, f.get_axes())
    leg = amp_ax.legend(loc='lower center', ncol=len(ranges) // 2, fancybox=True, frameon=True, shadow=False, framealpha=None, handlelength=0.3, numpoints=4)
    leg.get_title().set_fontsize(plt.rcParams['legend.fontsize'])
    cbar_ax = f.add_axes([0.95, 0.15, 0.03, 0.7])
    cbar = f.colorbar(mappable=mappable, cax=cbar_ax, boundaries=bd, norm=nm)#colorbar encoding distance
    #tick position is at the midpoint of ranges
    tick_pos = (bd[1:] + bd[:-1])/2
    cbar.set_ticks(tick_pos)
    cbar.set_ticklabels(["{:.2f}".format(r) for r in ranges])
    cbar.set_clim(ranges[0], ranges[-1])
    cbar.set_label('Distance from radar [m]')
    f.subplots_adjust(top=1)
    f.savefig(outputs[0])
    plt.close(f)


plot_figure_2(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)
