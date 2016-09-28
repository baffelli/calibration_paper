import matplotlib.pyplot as plt
import numpy as np
import pyrat.core.matrices as mat
import pyrat.fileutils.gpri_files as gpf
import pyrat.geo.geofun as geo
import pyrat.visualization.visfun as vf

def plot_reflectors(inputs, outputs, threads, config, params, wildcards):
    # load dem par
    dem_par = gpf.par_to_dict(inputs['dem_seg_par'])
    cov_par = gpf.par_to_dict(params['cov_par'])
    # T_mat = mat.coherencyMatrix(params['C_root'], params['C_root'] + '.par', gamma=True, bistatic=True,
    #                             basis='lexicographic').to_monostatic().lexicographic_to_pauli()
    # load coherency matrix
    C_mat_gc = np.array([gpf.load_binary(C + '_gc', dem_par['width']) for C in sorted(inputs['C'])]).transpose(
        (2, 1, 0))
    # load LUT
    lut = gpf.load_binary(inputs['lut'], dem_par['width'])
    #Compute validity mask
    invalid_mask = geo.gc_map_mask((cov_par['range_samples'], cov_par['azimuth_lines']), lut).T
    #Convert to right shape
    C_mat_gc = C_mat_gc.reshape((C_mat_gc.shape[0], C_mat_gc.shape[1], 4, 4))
    C_mat_gc = mat.coherencyMatrix(C_mat_gc, basis='lexicographic',
                                   bistatic=True).to_monostatic().lexicographic_to_pauli()
    T_mat_gc = C_mat_gc.boxcar_filter([5, 3])
    H_gc, A_gc, alpha_gc, beta_gc, p, w = T_mat_gc.cloude_pottier()
    # Compute pauli rgb
    rgb = vf.mask_zeros(T_mat_gc.pauli_image(k=0.3, sf=0.0005))

    with plt.style.context(config['style']):
        w, h = plt.rcParams['figure.figsize']
        alpha_fig = plt.figure(figsize=(w, w),)
        alpha_ax = alpha_fig.add_axes([0,0,1,1])
        H_fig = plt.figure(figsize=(w, w))
        H_ax = H_fig.add_axes([0,0,1,1])
        rgb_fig = plt.figure(figsize=(w, w))
        rgb_ax = rgb_fig.add_axes([0,0,1,1])
        # plot alpha
        alpha_mappable = alpha_ax.imshow(np.ma.array(np.rad2deg(alpha_gc)), mask=invalid_mask, vmin=0, vmax=90,
                                         interpolation='none')
        alpha_fig.colorbar(alpha_mappable, label=r'$\alpha [deg]$', orientation='horizontal', fraction=0.05, shrink=0.25, pad=0.02)
        alpha_ax.axis('off')
        kwdict = {'left':0, 'right':1, 'top':1, 'bottom':0, 'hspace':0, 'wspace':1}
        alpha_fig.subplots_adjust(**kwdict)
        # plot h
        H_mappable = H_ax.imshow(np.ma.array(H_gc), mask=invalid_mask, vmin=0, vmax=1, interpolation='none', cmap='gray')
        H_fig.colorbar(H_mappable, label=r'$H$', orientation='horizontal', fraction=0.05, shrink=0.35, pad=0.02)
        H_ax.axis('off')
        H_fig.subplots_adjust(**kwdict)
        # plot rgb paulo
        rgb_ax.imshow(rgb, interpolation='none')
        rgb_ax.axis('off')
        rgb_fig.subplots_adjust(**kwdict)
        H_fig.savefig(outputs['H'],  pad_inches=0)
        alpha_fig.savefig(outputs['alpha'], bbox_inches=0, pad_inches=0)
        rgb_fig.savefig(outputs['pauli'], bbox_inches=0, pad_inches=0)


plot_reflectors(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
                snakemake.wildcards)
