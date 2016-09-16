import pyrat.core.matrices as mat
import pyrat.core.polfun as pf
import matplotlib.pyplot as plt
import pyrat.visualization.visfun as vf
import pyrat.fileutils.gpri_files as gpf
import pyrat.geo.geofun as geo
import numpy as np
import osgeo.gdal as gdal
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import pyrat.core.corefun as cf
import json


def plot_reflectors(inputs, outputs, threads, config, params, wildcards):
    # load dem par
    dem_par = gpf.par_to_dict(inputs['dem_seg_par'])
    cov_par = gpf.par_to_dict(params['cov_par'])
    # load coherency matrix
    C_mat_gc = np.array([gpf.load_binary(C + '_gc', dem_par['width']) for C in sorted(inputs['C'])]).transpose(
        (1, 2, 0))
    print(C_mat_gc.shape)
    C_mat_gc = C_mat_gc.reshape((C_mat_gc.shape[0], C_mat_gc.shape[1], 4, 4))
    C_mat_gc = mat.coherencyMatrix(C_mat_gc, basis='lexicographic',
                                   bistatic=True).to_monostatic().lexicographic_to_pauli()
    T_mat_gc = C_mat_gc.boxcar_filter([5, 2])
    H_gc, A_gc, alpha_gc, beta_gc, p, w = T_mat_gc.cloude_pottier()
    invalid_mask = np.logical_or(T_mat_gc.span() == 0, T_mat_gc.span() == np.nan)
    # H_gc[T_mat_gc.span() < 1e-3] = np.nan
    # alpha_gc[T_mat_gc.span() < 1e-3] = np.nan
    # Compute pauli rgb
    rgb = vf.mask_zeros(C_mat_gc.pauli_image(k=0.5, sf=2))

    with plt.style.context(config['style']):
        w, h = plt.rcParams['figure.figsize']
        alpha_fig, alpha_ax = plt.subplots(figsize=(w, w))
        H_fig, H_ax = plt.subplots(figsize=(w, w))
        rgb_fig, rgb_ax = plt.subplots(figsize=(w, w))
        # plot alpha
        alpha_mappable = alpha_ax.imshow(np.rad2deg(alpha_gc), vmin=0, vmax=90, interpolation='none')
        alpha_fig.colorbar(alpha_mappable, label=r'$\alpha [deg]$', orientation='horizontal', shrink=0.5)
        alpha_ax.axis('off')
        alpha_fig.tight_layout()
        # plot h
        H_ax.imshow(H_gc, vmin=0, vmax=1, interpolation='none', cmap='gray')
        H_ax.axis('off')
        H_fig.tight_layout()
        # plot rgb paulo
        rgb_ax.imshow(rgb, interpolation='none')
        rgb_ax.axis('off')
        rgb_fig.tight_layout()
        H_fig.savefig(outputs['H'])
        alpha_fig.savefig(outputs['alpha'])
        rgb_fig.savefig(outputs['pauli'])


plot_reflectors(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
                snakemake.wildcards)
