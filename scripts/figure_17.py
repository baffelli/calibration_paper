

import pyrat.fileutils.gpri_files as gpf
import pyrat.core.matrices as mat
import pyrat.core.polfun as pf
import pyrat.core.corefun as cf
import matplotlib.pyplot as plt
import pyrat.gpri_utils.processors as proc
import pyrat.visualization.visfun as vf
from mpl_toolkits.mplot3d import Axes3D
import pyrat.gpri_utils.calibration as cal
import numpy as np
import csv


def figure_17(inputs, outputs, threads, config, params, wildcards):
    n_spec = 4#split the spectrum into 5 parts
    sw =(16,40)#search window
    # #list of reflectors
    refl_list = [ref for ref in params['ref'] if ref['type'] == "cubic" or ref['type'] == 'triangular']
    refl_list = sorted(refl_list, key=lambda x: x['ridx'])
    #Load matrix and convert to coherency
    slc_VV  = gpf.gammaDataset(inputs['slc_VV_par'], inputs['slc_VV'])
    slc_VV_split = proc.split_range_bandwidth(slc_VV, n_spec)
    f, (ax1,ax2) = plt.subplots(2,1, sharex=True, sharey=True)
    rgb1, *rest = vf.dismph(slc_VV_split[0])
    rgb2, *rest = vf.dismph(slc_VV_split[-1])
    ax1.imshow(rgb1)
    ax2.imshow(rgb2)
    plt.show()
    for ref in refl_list:
            ridx = ref['ridx']
            azidx = ref['azidx'] / slc_VV.azimuth_looks
            sl = (slice(ridx-sw[0]//2,ridx+sw[0]//2), slice(azidx-sw[1]//2,azidx+sw[1]//2))
            for slc in slc_VV_split:
                slc = cf.complex_interp(slc, (n_spec,1))
                ph_shift_VV, res_VV, r_sl_VV, meas_ph_VV, sim_ph_VV, rvp_VV, rcm_VV =  cal.measure_phase_center_location(slc_VV, ridx, azidx, sw=sw, aw=60)
                print(ph_shift_VV)
                plt.imshow(vf.dismph(slc[sl])[0])
                plt.show()


figure_17(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params, snakemake.wildcards)