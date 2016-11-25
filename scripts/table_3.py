


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

def tcr_recap(inputs, outputs, threads, config, params, wildcards):
    #list of reflectors
    refl_list = params['ref']
    #remove reflector used for calibration
    refl_list.pop(config['calibration']['reflector_index'])
    #Load matrix and convert to coherency
    C_root, ext = path.splitext(inputs['C_par'])
    C = mat.coherencyMatrix(C_root,  inputs['C_par'], gamma=True, bistatic=True, basis='lexicographic').boxcar_filter([3,3])
    # refl_list = [ref for ref in params['ref'] if ref['type'] == "cubic" or ref['type'] == 'triangular']
    refl_list = sorted(refl_list, key=lambda x: x['ridx'])
    with open(outputs[0], 'w+') as of:
        tabwrite = csv.writer(of, delimiter=',')
        tabwrite.writerow(
            ['name', 'slant range', 'HH-VV amplitude imbalance', 'HH-VV phase imbalance', 'Polarization purity', 'RCS', 'range_index', "azimuth_index"])
        f_arr = []
        HHVV_arr = []
        for ref in refl_list:
            idx_r = ref['ridx']
            idx_az = ref['azidx'] / C.GPRI_decimation_factor
            ptarg_zoom_C, rplot_C, azplot_C, mx_pos, resolution_dict, azvec, rvec = cf.ptarg(C[:,:,:,:], float(idx_r), float(idx_az), azwin=10, rwin=10, sw=(2,4))
            ptarg_zoom_C = mat.coherencyMatrix(ptarg_zoom_C, basis='lexicographic', bistatic=True)
            #Now estimate calibration parameters
            C_zoom = ptarg_zoom_C[mx_pos[0], mx_pos[1], :, :]
            r_sl = C.r_vec[idx_r]
            HHVV_phase = np.rad2deg(np.angle(C_zoom[0,3]))
            f = (C_zoom[3,3].real / C_zoom[0,0].real)**(1/4.0)
            purity = cf.dB(C_zoom[(0,0)].real) - cf.dB(C_zoom[(1,1)].real)
            C_sigma = C_zoom[0,0]
            RCS = C_sigma.real / cal.cr_rcs(ref['side'], C.radar_frequency, type=ref['type'])
            str = "{HHVV_phase},{f}".format(HHVV_phase=HHVV_phase, f=f)
            row = [ref['name'], r_sl, f, HHVV_phase, purity, RCS,  idx_r, idx_az]
            tabwrite.writerow(row)


tcr_recap(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params, snakemake.wildcards)