


import pyrat.core.matrices as mat
import pyrat.core.polfun as pf
import pyrat.core.corefun as cf
import matplotlib.pyplot as plt
import pyrat.visualization.visfun as vf
from mpl_toolkits.mplot3d import Axes3D
import pyrat.gpri_utils.calibration as cal
import numpy as np
import csv
import os.path as path
import pandas as pd

def tcr_recap(inputs, outputs, threads, config, params, wildcards):
    #list of reflectors
    with open(outputs[0], 'w+') as of:
        tabwrite = csv.writer(of, delimiter=',')
        tabwrite.writerow(
            ['RMS phase imbalance', 'RMS amplitude imbalance', 'RMS RCS bias'])
        residuals = pd.read_csv(inputs.res)
        rms_ph = np.mean(np.sqrt(residuals['HH-VV phase imbalance']**2))
        rms_f = np.mean(np.sqrt(residuals['HH-VV amplitude imbalance'] ** 2))
        tabwrite.writerow([rms_ph,rms_f, 0])
        # f_arr = []
        # HHVV_arr = []
        # for idx_r, idx_az, ty in refl_list_dec:
        #     ptarg_zoom_C, rplot_C, azplot_C, mx_pos, resolution_dict, azvec, rvec = cf.ptarg(C[:,:,:,:], float(idx_r), float(idx_az), azwin=10, rwin=10, sw=(2,4))
        #     ptarg_zoom_C = mat.coherencyMatrix(ptarg_zoom_C, basis='lexicographic', bistatic=True)
        #     #Now estimate calibration parameters
        #     C_zoom = ptarg_zoom_C[mx_pos[0], mx_pos[1], :, :]
        #     r_sl = C.r_vec[idx_r]
        #     HHVV_phase = np.rad2deg(np.angle(C_zoom[0,3]))
        #     f = (C_zoom[3,3].real / C_zoom[0,0].real)**(1/4.0)
        #     purity = cf.dB(C_zoom[(0,0)].real) - cf.dB(C_zoom[(1,1)].real)
        #     RCS = cf.dB(C_zoom[(0,0)].real)
        #     str = "{HHVV_phase},{f}".format(HHVV_phase=HHVV_phase, f=f)
        #     row = [r_sl, f, HHVV_phase, purity, RCS,  idx_r, idx_az]
        #     if ty == "t":
        #         f_arr.append(f)
        #         HHVV_arr.append(HHVV_phase)
        #     tabwrite.writerow(row)
        # # tabwrite.writerow(["RMS f", "RMS HHVV phase"])
        # # tabwrite.writerow([np.sqrt(np.mean(f**2)), np.sqrt(np.mean(HHVV_phase**2))])

tcr_recap(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params, snakemake.wildcards)