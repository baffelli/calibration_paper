import matplotlib.pyplot as plt
import pyrat.fileutils.gpri_files as gpf
import numpy as _np
import json
import pyrat.visualization.visfun as vf
import csv

import pyrat.gpri_utils.calibration as cal

def plot_azimuth_phase(inputs, outputs, threads, config, params):
    #import the slcs
    slc = gpf.gammaDataset(inputs['slc'] + '.par', inputs['slc'])
    #load list of reflectors
    refl_list = config['list_of_reflectors']
    print(refl_list)
    #Create figure
    with plt.style.context(config['style']):
        f, (phase_ax, amp_ax) = plt.subplots(2, sharex=True)
        sorted_by_range = sorted(refl_list, key=lambda tup: tup[0])
        sorted_by_range = [ref for ref in sorted_by_range if ref[-1] == "t"]
        for ridx, azidx, *rest in sorted_by_range:
                #Slice the slc
                slc_sl = (ridx, slice(azidx - params['ws'] / 2, azidx + params['ws']/2))
                subimage = slc[slc_sl]
                #Determine true maximum in the slice
                max_idx = _np.argmax(_np.abs(subimage))
                #Determine the shift
                shift = subimage.shape[0] / 2 - max_idx
                #Extract a new slice that has been shifted
                slc_slc_new = slc_sl = (ridx, slice(azidx - shift - params['ws'] / 2, azidx - shift + params['ws']/2))
                reflector_slice = slc[slc_slc_new]
                #Slice slc
                #Azimuth angle vector for plot
                az_vec = slc.GPRI_az_angle_step[0] * _np.arange(-len(reflector_slice)/2
                                                                                  ,len(reflector_slice)/2)
                refl_ph = _np.angle(reflector_slice)
                refl_ph = _np.unwrap(refl_ph)
                refl_ph -= refl_ph[reflector_slice.shape[0]/2]
                max_phase = refl_ph[reflector_slice.shape[0]/2]
                refl_amp = (_np.abs(reflector_slice))**2
                r_sl = slc.r_vec[ridx]
                line, = phase_ax.plot(az_vec,_np.rad2deg(refl_ph), label=r"r={} m".format(round(r_sl)))
                amp_ax.plot(az_vec, refl_amp/refl_amp[reflector_slice.shape[0]/2])
                #Plot line for beamwidth
                phase_ax.set_ylim(-30,30)
                phase_ax.axvline(0.2, color='red', ls='--')
                phase_ax.axvline(-0.2, color='red', ls='--')
                phase_ax.yaxis.set_label_text(r'Phase [deg]')
                # phase_ax.xaxis.set_label_text(r'azimuth angle from maximum [deg]')
                amp_ax.axvline(0.2, color='red', ls='--')
                amp_ax.axvline(-0.2, color='red', ls='--')
                amp_ax.yaxis.set_label_text(r'Relative Intensity')
                amp_ax.xaxis.set_label_text(r'azimuth angle from maximum [deg]')
                phase_ax.legend()
                f.savefig(outputs['plot'])
                plt.close(f)
plot_azimuth_phase(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params)