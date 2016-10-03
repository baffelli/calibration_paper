import matplotlib.pyplot as plt
import numpy as np
import pyrat.fileutils.gpri_files as gpf
import pyrat.core.corefun as cf

def plot_figure_2(inputs, outputs, threads, config, params, wildcards):
    plt.style.use(inputs['style'])
    slc = gpf.gammaDataset(inputs['VV'] + '.par', inputs['VV'])
    f, (phase_ax, amp_ax) = plt.subplots(2, sharex=True)
    sorted_by_range = sorted(params['reflectors'], key=lambda tup: tup[0])
    sorted_by_range = [ref for ref in sorted_by_range if ref[-1] == "t"]
    for ridx, azidx, *rest in sorted_by_range:
        #Extract maximum around the supposed position
        ridx_max, az_idx_max = cf.maximum_around(np.abs(slc), [ridx, azidx], [2,5])
        # Slice the slc
        slc_sl = (ridx_max, slice(azidx - params['ws'] / 2, azidx + params['ws'] / 2))
        reflector_slice = slc[slc_sl]
        max_slc = slc[ridx_max, az_idx_max]
        # Azimuth angle vector for plot
        az_vec = slc.GPRI_az_angle_step[0] * np.arange(-len(reflector_slice) / 2
                                                        , len(reflector_slice) / 2)
        refl_ph = np.angle(reflector_slice) - np.angle(max_slc)
        refl_ph = np.unwrap(refl_ph)
        refl_amp = (np.abs(reflector_slice)) ** 2
        r_sl = slc.r_vec[ridx]
        # line, = phase_ax.plot(az_vec, _np.rad2deg(refl_ph), label=r"r={} m".format(round(r_sl)))
        amp_ax.plot(az_vec, refl_amp / refl_amp[reflector_slice.shape[0] / 2], label=r"r={} m".format(round(r_sl)))
        # Plot line for beamwidth
        phase_ax.set_ylim(-30, 30)
        phase_ax.axvline(0.2, color='red', ls='--')
        phase_ax.axvline(-0.2, color='red', ls='--')
        phase_ax.yaxis.set_label_text(r'Phase [deg]')
        # phase_ax.xaxis.set_label_text(r'azimuth angle from maximum [deg]')
        amp_ax.axvline(0.2, color='red', ls='--')
        amp_ax.axvline(-0.2, color='red', ls='--')
        amp_ax.yaxis.set_label_text(r'Relative Intensity')
        amp_ax.xaxis.set_label_text(r'azimuth angle from maximum [deg]')
        amp_ax.legend()
        f.savefig(outputs[0])
        plt.close(f)

plot_figure_2(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)