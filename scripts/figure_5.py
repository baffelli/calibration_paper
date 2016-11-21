import os.path as path

import matplotlib.pyplot as plt
import numpy as np
import pyrat.core.corefun as cf
import pyrat.core.matrices as mat
import pyrat.core.polfun as pf

from matplotlib.ticker import MultipleLocator

from mpl_toolkits.mplot3d import Axes3D

# Return the decimated azimuth position
def az_idx(ds, idx):
    print(ds.azimuth_looks)
    return idx / ds.GPRI_decimation_factor


def plot_figure_5(inputs, outputs, threads, config, params, wildcards):
    plt.style.use(inputs['style'])
    fig_w, fig_h = plt.rcParams['figure.figsize']
    f, fig_arr = plt.subplots(2, 2, subplot_kw=dict(projection='3d'), figsize=(fig_w * 2, fig_h * 2))#two column figure

    C_par = inputs["C_par"]
    C_root, ext = path.splitext(C_par)
    C_cal_par = inputs["C_cal_par"]
    C_cal_root, ext = path.splitext(C_cal_par)
    C = mat.coherencyMatrix(C_root, C_par, gamma=True, bistatic=True,
                            basis='lexicographic').to_monostatic().boxcar_filter([3, 3])
    C_cal = mat.coherencyMatrix(C_cal_root, C_cal_par, gamma=True, bistatic=True,
                                basis='lexicographic').to_monostatic().boxcar_filter([3, 3])

    #find the maxmium
    ridx = params.ref['ridx']
    azidx = az_idx(C, params.ref['azidx'])
    mx_idx = cf.maximum_around(np.abs(C), [ridx,azidx, 0,0], [2,4,1,1])
    #compute signature
    co_sig_cal, x_sig_cal, psi, chi = pf.pol_signature(C_cal[mx_idx[0],mx_idx[1]], n_points=300)
    co_sig, x_sig, psi, chi = pf.pol_signature(C[mx_idx[0], mx_idx[1]],
                                                       n_points=300)
    #extract the maximum on the
    pp, cc = (np.rad2deg(psi), np.rad2deg(chi))
    # plot uncalibrated
    fig_arr[0, 0].plot_surface(pp, cc, co_sig, cmap='RdBu_r', lw=0.2)
    fig_arr[0, 1].plot_surface(pp, cc, x_sig, cmap='RdBu_r', lw=0.2)
    #plot calibrated
    fig_arr[1, 0].plot_surface(pp, cc, co_sig_cal, cmap='RdBu_r', lw=0.2)
    fig_arr[1, 1].plot_surface(pp, cc, x_sig_cal, cmap='RdBu_r', lw=0.2)
        # fig_arr[ref_count, 0].auto_scale_xyz([-90, 90], [-45, 45], [0, 1])
    #Set label etc https://dawes.wordpress.com/2014/06/27/publication-ready-3d-figures-from-matplotlib/
    for row in fig_arr:
        for ax in row:
            ax.view_init(elev=45, azim=45)
            ax.xaxis.set_major_locator(MultipleLocator(45))
            ax.yaxis.set_major_locator(MultipleLocator(45))
            ax.zaxis.set_major_locator(MultipleLocator(0.5))
            ax.set_ylabel(r' $\chi$ [deg]', labelpad=0)
            ax.set_xlabel(r' $\psi$ [deg]',  labelpad=0)
            ax.set_zlabel(r'Power',  labelpad=0)
            ax.dist = 12.5
            ax.auto_scale_xyz([-90, 90], [-45, 45], [0, 1])
    f.subplots_adjust(wspace=0.1, hspace=0.1, bottom=0.25, top=1, left=0.2)
    f.suptitle(params.ref['name'])
    f.savefig(outputs[0])


plot_figure_5(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)
