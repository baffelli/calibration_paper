


import pyrat.core.matrices as mat
import pyrat.core.polfun as pf
import pyrat.core.corefun as cf
import matplotlib.pyplot as plt
import pyrat.visualization.visfun as vf
from mpl_toolkits.mplot3d import Axes3D
import pyrat.gpri_utils.calibration as cal
import numpy as np
import csv

def tcr_recap(inputs, outputs, threads, config, params, wildcards):
    #list of reflectors
    refl_list = config['list_of_reflectors']
    refl_list_dec = [ [ridx, azidx / config['range_compression']['dec']] for ridx, azidx in refl_list]
    #Load matrix and convert to coherency
    C = mat.coherencyMatrix(params['C_root'],  params['C_root'] + '.par', gamma=True, bistatic=True, basis='lexicographic').boxcar_filter([3,3])
    with open(outputs['cal_report'], 'w+') as of:
        tabwrite = csv.writer(of, delimiter=',')
        tabwrite.writerow(
            ['HH-VV phase imbalance', 'HH-VV amplitude imbalance', 'Polarisation purity', 'RCS', 'slant range'])
        for idx_r, idx_az in refl_list_dec:
            ptarg_zoom_C, rplot_C, azplot_C, mx_pos, resolution_dict = cf.ptarg(C[:,:,:,:], float(idx_r), float(idx_az), azwin=10, rwin=10, sw=2)
            ptarg_zoom_C = mat.coherencyMatrix(ptarg_zoom_C, basis='lexicographic', bistatic=True)
            #Now estimate calibration parameters
            r_sl = C.r_vec[idx_r]
            C = ptarg_zoom_C[mx_pos]
            HHVV_phase = np.rad2deg(np.angle(C[0,3]))
            f = (C[3,3].real / C[0,0].real)**(1/4.0)
            purity = cf.dB(C[mx_pos + (0,0)].real) - cf.dB(C[mx_pos + (1,1)].real)
            RCS = cf.dB(C[mx_pos + (0,0)].real)
            str = "{HHVV_phase},{f}".format(HHVV_phase=HHVV_phase, f=f)
            row = [HHVV_phase, f, purity, RCS, r_sl]
            print(row)
            tabwrite.writerow(row)

tcr_recap(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params, snakemake.wildcards)