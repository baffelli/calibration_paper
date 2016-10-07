

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

def table_2(inputs, outputs, threads, config, params, wildcards):
    #list of reflectors
    refl_list = [ref for ref in params['ref'] if ref['type'] == "cubic" or ref['type'] == 'triangular']
    refl_list = sorted(refl_list, key=lambda x: x['ridx'])
    #Load matrix and convert to coherency
    slc = gpf.gammaDataset(inputs['slc_par'], inputs['slc'])
    #list of results
    res_list = []
    with open(outputs[0], 'w+') as of:
        tabwrite = csv.writer(of, delimiter=',')
        tabwrite.writerow(
            ["rph", 'res', 'rsl', "ridx", "azidx"])
        for ref in refl_list:
            ridx = ref['ridx']
            azidx = ref['azidx'] / slc.azimuth_looks
            ph_shift, res, r_sl, meas_ph, sim_ph =  cal.measure_phase_center_location(slc, ridx, azidx, sw=(2,25), aw=60)
            row = [ph_shift, res, r_sl, ridx, azidx]
            res_list.append(row)
            print(row)
            tabwrite.writerow(row)


table_2(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params, snakemake.wildcards)