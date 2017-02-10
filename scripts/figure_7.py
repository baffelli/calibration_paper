import os.path as path

import matplotlib.pyplot as plt
import numpy as np
import pyrat.core.corefun as cf
import pyrat.core.matrices as mat

import pyrat.visualization.visfun as vf

import pyrat.geo.geofun as geo

import cartopy.crs as ccrs
import cartopy

import osr

import gdal

import pyrat.fileutils.gpri_files as gpf

import collections
from scipy.special import expit

def annotate_axis(ax, params):
    # box = dict(boxstyle='round', fc='white', alpha=1)
    pos_list = collections.defaultdict(lambda:(30,0))
    pos_list.update( {'CR1': (30, 10), 'CR3':(20,0), 'CR4':(30,-10),'CR6': (10,-20)})
    for ref in params['ref']:#iterate reflector and positions
        ref_pos = ref['geo_coord']
        ax.plot(ref_pos[0], ref_pos[1], mec=ref['marker_color'], marker='o', mfc='none', ms=10)
        ax.annotate(
            ref['name'],
            xy = ref_pos, xytext = pos_list[ref['name']],
            textcoords = 'offset points', ha = 'right', va = 'bottom', color='white')


def plot_figure_7(inputs, outputs, threads, config, params, wildcards):
    #load data
    C_cal_par = inputs["C_cal_par"]
    C_cal_root, ext = path.splitext(C_cal_par)
    C_cal = mat.coherencyMatrix(C_cal_root, C_cal_par, gamma=True, bistatic=True,
                                basis='lexicographic').lexicographic_to_pauli()

    #Load shadow map
    sh_map = gpf.gammaDataset(inputs.dem_seg_par, inputs.sh_map, dtype=gpf.type_mapping['UCHAR'])

    #load map
    map_ds = gdal.Open(inputs.map)
    #Load dem segment
    dem_ds = gdal.Open(inputs.dem_seg)
    #Create spatial reference
    ref = osr.SpatialReference()
    ref.ImportFromWkt(map_ds.GetProjection())

    #create figure
    plt.style.use(inputs['style'])
    fig_w, fig_h = plt.rcParams['figure.figsize']
    f, ax = plt.subplots(figsize=(fig_w, fig_h*2))

    #Geocode
    LUT = geo.GeocodingTable(inputs.dem_seg_par, inputs.LUT, inputs.C_cal_par, inputs.LUT_inv)

    #Decimate reflector location
    for ref in params['ref']:
        ref['azidx_dec'] = C_cal.azidx_dec(ref['azidx'])
        ref['geo_coord'] = LUT.radar_coord_to_geo_coord([ref['ridx'], ref['azidx_dec']])[0]
        ref['marker_color'] = 'orange' if ref['name'] == params['ref'][config['calibration']['reflector_index']]['name'] else '#43a2ca'
    # #Geocode
    C_cal_geo = LUT.geocode_data(C_cal)
    C_cal_geo = C_cal_geo.boxcar_filter([3,3])
    C_cal_rgb = C_cal_geo.pauli_image(k=0.4, sf=0.25, common=False, peak=False)
    #Save RGB as a tif
    C_cal_alpha = np.ones(C_cal_rgb[:,:,0].shape)
    #Mask Nans
    C_cal_span = C_cal_geo.span()
    C_cal_alpha = np.ones_like(C_cal_span)
    # C_cal_alpha[(C_cal_geo.span()/np.nanmax(C_cal_span)) < 1e-6] = 0
    C_cal_alpha[np.isnan(C_cal_geo[:,:,0,0])] = 0
    #Convert to rgba
    C_cal_rgb = np.dstack((C_cal_rgb,C_cal_alpha))
    #Save as rif
    geo.raster_to_geotiff(C_cal_geo.pauli_image(k=0.4, sf=0.25, common=False).transpose((1,0,2)), dem_ds, outputs.tif)
    #Plot map
    ext1=geo.get_ds_extent(map_ds)
    ext2=LUT.get_extent()
    ax.imshow(map_ds.ReadAsArray().transpose((1, 2, 0)), extent=ext1)
    ax.imshow(C_cal_rgb.transpose(1,0,2), extent=ext2, alpha=0.95)
    #annotate scatters
    annotate_axis(ax, params)
    ext = LUT.get_geocoded_extent(C_cal)
    ax.set_xlim([606637,609660])
    ax.set_xlim(ext[0:2])
    ax.set_ylim(ext[2::])
    ax.set_aspect(1)
    # f.subplots_adjust(left=0,right=1)
    ax.axis('off')
    plt.show()
    f.savefig(outputs['paper_fig'], bbox_inches='tight', pad_inches=0)
    #Add axis
    f1, ax1 =  plt.subplots(figsize=( fig_w, fig_w))
    ax1.imshow(C_cal_rgb.transpose(1, 0, 2), extent=ext2, aspect=1)
    ax1.axis('off')
    # annotate scatters
    annotate_axis(ax1, params)
    f1.subplots_adjust(left=0, right=1)
    f1.savefig(outputs['pres_fig'])



plot_figure_7(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)
