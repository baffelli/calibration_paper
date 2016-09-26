


import pyrat.core.matrices as mat
import pyrat.core.polfun as pf
import pyrat.core.corefun as cf
import matplotlib.pyplot as plt
import pyrat.visualization.visfun as vf
from mpl_toolkits.mplot3d import Axes3D
import pyrat.gpri_utils.calibration as cal
import numpy as np
import csv
def est_residuals(inputs, outputs, threads, config, params, wildcards):
    with plt.style.context(config['style']):
        C = mat.coherencyMatrix(params['C_root'],  params['C_root'] + '.par', gamma=True, bistatic=True, basis='lexicographic').boxcar_filter([3,3])
        #Find the maxium by interpolation
        #Compute the ptarg response for the C matrix
        ptarg_zoom_C, rplot_C, azplot_C, mx_pos, resolution_dict = cf.ptarg(C[:,:,:,:], float(wildcards['ridx']), float(wildcards['azidx']), azwin=20, rwin=10, sw=(2,4))
        ptarg_zoom_C = mat.coherencyMatrix(ptarg_zoom_C, basis='lexicographic', bistatic=True)
        C_tri = ptarg_zoom_C[mx_pos]
        #Now estimate calibration parameters
        r_sl = C.r_vec[int(wildcards['ridx'])]
        HHVV_phase = np.rad2deg(np.angle(C_tri[0,3]))
        f = (C_tri[3,3].real / C_tri[0,0].real)**(1/4.0)
        purity = cf.dB(ptarg_zoom_C[mx_pos + (0,0)].real) - cf.dB( ptarg_zoom_C[mx_pos + (1,1)].real)
        RCS = cf.dB(ptarg_zoom_C[mx_pos + (0,0)].real)
        str = "{HHVV_phase},{f}".format(HHVV_phase=HHVV_phase, f=f)
        with open(outputs['est_parameters'],'w+') as of:
            tabwrite = csv.writer(of, delimiter=',')
            tabwrite.writerow(['HH-VV phase imbalance', 'HH-VV amplitude imbalance', 'Polarisation purity', 'RCS', 'slant range'])
            tabwrite.writerow([HHVV_phase, f, purity, RCS, r_sl])
est_residuals(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params, snakemake.wildcards)