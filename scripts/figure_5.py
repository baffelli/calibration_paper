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
    return idx / ds.azimuth_looks


def plot_figure_5(inputs, outputs, threads, config, params, wildcards):
    plt.style.use(inputs['style'])
    fig_w, fig_h = plt.rcParams['figure.figsize']
    f, fig_arr = plt.subplots(4, 2, subplot_kw=dict(projection='3d'), figsize=(fig_w, fig_h * 4))

    C_par = inputs["C_par"]
    C_root, ext = path.splitext(C_par)
    C_cal_par = inputs["C_cal_par"]
    C_cal_root, ext = path.splitext(C_cal_par)
    C = mat.coherencyMatrix(C_root, C_par, gamma=True, bistatic=True,
                            basis='lexicographic').to_monostatic().boxcar_filter([3, 3])
    C_cal = mat.coherencyMatrix(C_cal_root, C_cal_par, gamma=True, bistatic=True,
                                basis='lexicographic').to_monostatic().boxcar_filter([3, 3])
    for ref_count, ref_idx in zip([0,2],[1, 5]):
        #find the maxmium
        ridx = params.ref[ref_idx][0]
        azidx = az_idx(C, params.ref[ref_idx][1])

        mx_idx = cf.maximum_around(np.abs(C), [ridx,azidx, 0,0], [2,4,1,1])
        #compute signature
        co_sig_cal, x_sig_cal, psi, chi = pf.pol_signature(C_cal[mx_idx[0],mx_idx[1]], n_points=300)
        co_sig, x_sig, psi, chi = pf.pol_signature(C[mx_idx[0], mx_idx[1]],
                                                           n_points=300)
        #extract the maximum on the
        pp, cc = (np.rad2deg(psi), np.rad2deg(chi))
        # plot copol
        fig_arr[ref_count, 0].plot_surface(pp, cc, co_sig, cmap='RdBu_r', lw=0.2)
        fig_arr[ref_count, 1].plot_surface(pp, cc, x_sig, cmap='RdBu_r', lw=0.2)
        fig_arr[ref_count  + 1, 0].plot_surface(pp, cc, co_sig_cal, cmap='RdBu_r', lw=0.2)
        fig_arr[ref_count  + 1, 1].plot_surface(pp, cc, x_sig_cal, cmap='RdBu_r', lw=0.2)
        fig_arr[ref_count, 0].auto_scale_xyz([-90, 90], [-45, 45], [0, 1])
    #Set label etc https://dawes.wordpress.com/2014/06/27/publication-ready-3d-figures-from-matplotlib/
    for row in fig_arr:
        for ax in row:
            ax.view_init(elev=70, azim=30)
            ax.xaxis.set_major_locator(MultipleLocator(45))
            ax.yaxis.set_major_locator(MultipleLocator(45))
            ax.zaxis.set_major_locator(MultipleLocator(0.5))
            ax.set_ylabel(r' $\chi$ [deg]', linespacing=0.1)
            ax.set_xlabel(r' $\psi$ [deg]', linespacing=0.1)
            ax.set_zlabel(r'Power', linespacing=0.1)
            ax.xaxis.labelpad = 0.1
            ax.yaxis.labelpad = 0.1
            ax.zaxis.labelpad = 0,1
            # ax.xaxis.line.set_lw(0.5)
            # ax.yaxis.line.set_lw(0.5)
            # ax.zaxis.line.set_lw(0.5)
            # ax.distance = 13
            # ax.grid(False)
            # ax.xaxis.pane.set_edgecolor('black')
            # ax.yaxis.pane.set_edgecolor('black')
            # ax.xaxis.pane.fill = False
            # ax.yaxis.pane.fill = False
            # ax.zaxis.pane.fill = False
            ax.xaxis._axinfo['tick']['inward_factor'] = 0.4
            ax.xaxis._axinfo['tick']['outward_factor'] = 0
            ax.yaxis._axinfo['tick']['inward_factor'] = 0.4
            ax.yaxis._axinfo['tick']['outward_factor'] = 0
            ax.zaxis._axinfo['tick']['inward_factor'] = 0.4
            ax.zaxis._axinfo['tick']['outward_factor'] = 0
            # [t.set_va('center') for t in ax.get_yticklabels()]
            # [t.set_ha('left') for t in ax.get_yticklabels()]
            # [t.set_va('center') for t in ax.get_xticklabels()]
            # [t.set_ha('right') for t in ax.get_xticklabels()]
            # [t.set_va('center') for t in ax.get_zticklabels()]
            # [t.set_ha('left') for t in ax.get_zticklabels()]
            # print(ax.xaxis._axinfo)

    f.subplots_adjust(wspace=0.2, hspace=0.2, bottom=0.2)
    f.savefig(outputs[0])


plot_figure_5(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)
