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

from scipy.special import expit

# Return the decimated azimuth position
def az_idx(ds, idx):
    return idx / ds.azimuth_looks


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
    #Create spatial reference
    ref = osr.SpatialReference()
    ref.ImportFromWkt(map_ds.GetProjection())

    #create figure
    plt.style.use(inputs['style'])
    fig_w, fig_h = plt.rcParams['figure.figsize']
    f, ax = plt.subplots(figsize=(2 * fig_w,  2 * fig_w))

    #Geocode
    LUT = geo.GeocodingTable(inputs.dem_seg_par, inputs.LUT)

    #Decimate reflector location
    ref_dec = np.array([[ref['ridx'], ref['azidx'] // C_cal.GPRI_decimation_factor] for ref in params['ref']], dtype=int)
    ref_dec_geo = np.array([LUT.dem_coord_to_geo_coord(LUT.radar_coord_to_dem_coord(ref_pos)) for ref_pos in ref_dec])

    #Geocode
    C_cal_geo = LUT.geocode_data(C_cal)
    C_cal_geo = C_cal_geo.boxcar_filter([3,3])
    C_cal_rgb = C_cal_geo.pauli_image(k=0.4, sf=0.25, common=False, peak=False)
    C_cal_alpha = np.ones(C_cal_rgb[:,:,0].shape)
    #Mask Nans
    C_cal_alpha[np.isnan(C_cal_geo[:,:,0,0])] = 0
    #Convert to rgba
    C_cal_rgb = np.dstack((C_cal_rgb,C_cal_alpha))
    #Plot map
    ext1=geo.get_ds_extent(map_ds)
    ext2=LUT.get_extent()
    ax.imshow(map_ds.ReadAsArray().transpose((1, 2, 0)), extent=ext1)
    ax.imshow(C_cal_rgb.transpose(1,0,2), extent=ext2,)
    marker_color = ['orange' if i == config['calibration']['reflector_index'] else '#43a2ca' for i in range(len(ref_dec))]
    ax.scatter(ref_dec_geo[:,0], ref_dec_geo[:,1], edgecolors=marker_color, s=90, facecolors='none',
               linewidths=0.5)
    #annotate scatter
    pos_list = {'Chutzen': (45,20), 'Hindere Chlapf':(90,-10), 'Bifang':(10,10), 'Türle':(10,-20), 'Simmleremoos 1': (55,15),'Simmleremoos 2': (95,-20)}
    for ref, ref_pos in zip(params['ref'], ref_dec_geo):#iterate reflector and positions
        plt.annotate(
            ref['name'],
            xy = ref_pos, xytext = pos_list[ref['name']],
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0', color='grey'), color='white')
    ext = [606630,610072,188060,192957]
    ax.set_xlim(ext[0:2])
    ax.set_ylim(ext[2:])
    f.subplots_adjust(left=0,right=1)
    plt.show()

    # # C_cal_geo, x_vec, y_vec, lut, direct_lut = geo.geocode_image(C_cal, 2)
    # #Decimate reflector location
    #
    # #Compute locations in radar coordinates
    # ref_gc = vf.bilinear_interpolate(direct_lut.T, ref_dec[:,0],ref_dec[:,1])
    #
    # ax.imshow(C_cal_rgb.transpose([1,0,2])[::,::-1], aspect=1)
    # ax.scatter(C_cal_geo.shape[0] - np.imag(ref_gc), np.real(ref_gc), edgecolors=marker_color, s=25, facecolors='none', linewidths=0.5)
    ax.axis('off')
    f.savefig(outputs[0],dpi=1200)


plot_figure_7(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)
