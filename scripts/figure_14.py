import os.path as path

import matplotlib.pyplot as plt
import numpy as np
import pyrat.core.corefun as cf
import pyrat.core.matrices as mat

import pyrat.visualization.visfun as vf

import pyrat.core.polfun as pf

import pyrat.geo.geofun as geo

import cartopy.crs as ccrs
import cartopy

import osr

import gdal

import pyrat.fileutils.gpri_files as gpf

from scipy.special import expit

from matplotlib import gridspec



def plot_figure_14(inputs, outputs, threads, config, params, wildcards):
    #load data
    C_cal_par = inputs["C_cal_par"]
    C_cal_root, ext = path.splitext(inputs['C'][0])
    C_cal = mat.coherencyMatrix(C_cal_root, C_cal_par, gamma=True, bistatic=True,
                                basis='lexicographic').lexicographic_to_pauli().to_monostatic()

    #load map
    map_ds = gdal.Open(inputs.map)
    #Create spatial reference
    ref = osr.SpatialReference()
    ref.ImportFromWkt(map_ds.GetProjection())


    #Geocode
    LUT = geo.GeocodingTable(inputs.dem_seg_par, inputs.LUT)
    C_cal_geo = LUT.geocode_data(C_cal).boxcar_filter([5,2])
    #Compute decomposition
    H_geo, anisotropy, alpha_geo, beta_m_geo, p_geo, w_geo = C_cal_geo.cloude_pottier()
    #Geocode
    # alpha_geo  = LUT.geocode_data(alpha)
    # H_geo = LUT.geocode_data(H)
    #RGB composite
    plot_kw = {'k':0.2, 'sf':1.2, 'coherence':False, 'min_val':0, 'max_val':np.pi/2, 'black_background':False, 'type':'increasing'}
    span_geo = LUT.geocode_data(C_cal.span())
    rgb, norm, pal = vf.dismph(span_geo*np.exp(1j*alpha_geo), **plot_kw)
    pal, pal_ext = vf.dismph_palette(span_geo*np.exp(1j*alpha_geo), **plot_kw)
    #create figure
    plt.style.use(inputs['style'])
    fig_w, fig_h = plt.rcParams['figure.figsize']
    #Determine aspect and extent
    ext2=LUT.get_geocoded_extent(alpha_geo)
    ext2 = [ext2[-2], ext2[-1], ext2[0], ext2[1]]
    asp = vf.fixed_aspect(ext2, 1)
    #Add axis for H
    f_H = plt.figure(figsize=(2 * fig_w, 2 * fig_w))
    ax_H = f_H.add_subplot(111)
    ax_H.axis('off')
    map_H = ax_H.imshow(H_geo.T, cmap='gray')
    plt.colorbar(map_H)
    f_H.savefig(outputs['H_fig'])
    #Axis for alpha
    f_alpha = plt.figure(figsize=(2 * fig_w, 2 * fig_w))
    ax_alpha = f_alpha.add_subplot(111)
    ax_alpha.axis('off')
    map_alpha = ax_alpha.imshow(alpha_geo.T, extent=ext2, aspect=asp)
    plt.colorbar(map_alpha)
    f_alpha.savefig(outputs['alpha_fig'])
    #composite
    f1 =  plt.figure(figsize
                     =(2 * fig_w,  2 * fig_w))
    gs = gridspec.GridSpec(*(2, 3), height_ratios=[0.95, 0.05])
    ax1 = f1.add_subplot(gs[0, ::])
    cax = f1.add_subplot(gs[1,1])
    ax1.imshow(rgb.transpose(1, 0, 2), extent=ext2, aspect=asp)
    ax1.axis('off')
    cax.imshow(pal, extent=pal_ext, aspect=asp)
    cax.set_ylabel(r'$\alpha$')
    cax.set_xlabel(r'Intensity')
    cax.grid(b=False)
    cax.set_yticks([0, np.pi/4, np.pi/2])
    cax.set_xticks([0, pal_ext[1]/2, pal_ext[1]])
    cax.set_yticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$"])
    gs.update(hspace=0.05)
    f1.subplots_adjust(left=0, right=1)
    f1.savefig(outputs['rgb_fig'])


plot_figure_14(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)
