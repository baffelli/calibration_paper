


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
    refl_list_dec = sorted([ [ridx, azidx / config['range_compression']['dec'], ty] for ridx, azidx, ty in refl_list],key=lambda x: x[0])
    #Load matrix and convert to coherency
    C = mat.coherencyMatrix(params['C_root'],  params['C_root'] + '.par', gamma=True, bistatic=True, basis='lexicographic').boxcar_filter([3,3])
    with open(outputs['cal_report'], 'w+') as of:
        tabwrite = csv.writer(of, delimiter=',')
        tabwrite.writerow(
            ['HH-VV phase imbalance', 'HH-VV amplitude imbalance', 'Polarisation purity', 'RCS', 'slant range', 'range_index', "azimuth_index"])
        f_arr = []
        HHVV_arr = []
        for idx_r, idx_az, ty in refl_list_dec:
            ptarg_zoom_C, rplot_C, azplot_C, mx_pos, resolution_dict = cf.ptarg(C[:,:,:,:], float(idx_r), float(idx_az), azwin=10, rwin=10, sw=(2,4))
            ptarg_zoom_C = mat.coherencyMatrix(ptarg_zoom_C, basis='lexicographic', bistatic=True)
            #Now estimate calibration parameters
            C_zoom = ptarg_zoom_C[mx_pos[0], mx_pos[1], :, :]
            r_sl = C.r_vec[idx_r]
            HHVV_phase = np.rad2deg(np.angle(C_zoom[0,3]))
            f = (C_zoom[3,3].real / C_zoom[0,0].real)**(1/4.0)
            purity = cf.dB(C_zoom[(0,0)].real) - cf.dB(C_zoom[(1,1)].real)
            RCS = cf.dB(C_zoom[(0,0)].real)
            str = "{HHVV_phase},{f}".format(HHVV_phase=HHVV_phase, f=f)
            row = [HHVV_phase, f, purity, RCS, r_sl, idx_r, idx_az]
            if ty == "t":
                f_arr.append(f)
                HHVV_arr.append(HHVV_phase)
            print(row)
            tabwrite.writerow(row)
        tabwrite.writerow(["RMS f", "RMS HHVV phase"])
        tabwrite.writerow([np.sqrt(np.mean(f**2)), np.sqrt(np.mean(HHVV_phase**2))])

tcr_recap(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params, snakemake.wildcards)