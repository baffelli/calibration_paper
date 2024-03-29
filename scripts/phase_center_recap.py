

import pyrat.fileutils.gpri_files as gpf
import pyrat.core.matrices as mat
import pyrat.core.polfun as pf
import pyrat.core.corefun as cf
import matplotlib.pyplot as plt
import pyrat.visualization.visfun as vf
from mpl_toolkits.mplot3d import Axes3D
import pyrat.gpri_utils.calibration as cal
import numpy as np
import csv

def phc_recap(inputs, outputs, threads, config, params, wildcards):
    #list of reflectors
    refl_list = [ref for ref in config['list_of_reflectors'] if ref[-1]=="t"]
    refl_list = sorted([ [ridx, azidx, ty] for ridx, azidx, ty in refl_list],key=lambda x: x[0])
    #Load matrix and convert to coherency
    slc = gpf.gammaDataset(inputs['slc_par'], inputs['slc'])
    #list of results
    res_list = []
    with open(outputs['phase_report'], 'w+') as of:
        tabwrite = csv.writer(of, delimiter=',')
        tabwrite.writerow(
            ['phase_center_offset', 'residual', 'slant range', 'ridx', 'azidx'])
        for idx_r, idx_az, ty in refl_list:
            ph_shift, res, r_sl, meas_ph, sim_ph =  cal.measure_phase_center_location(slc, idx_r, idx_az, aw=20)
            row = [ph_shift, res, r_sl, idx_r, idx_az]
            res_list.append(row)
            print(row)
            tabwrite.writerow(row)
        tabwrite.writerow(['average r_ph', 'weighted average r_ph'])
        res_list = np.array(res_list)
        r_ph_opt = np.mean(res_list[:,0])
        r_ph_opt_w = np.average(res_list[:, 0], weights = 1/res_list[:,1])
        tabwrite.writerow([r_ph_opt, r_ph_opt_w])
        #Compute average
        #     ax.plot(meas_ph,'r')
        #     ax.plot(sim_ph,'g')
        # plt.show()

phc_recap(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params, snakemake.wildcards)