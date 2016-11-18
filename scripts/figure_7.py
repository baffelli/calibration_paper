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

# Return the decimated azimuth position
def az_idx(ds, idx):
    return idx / ds.azimuth_looks


def plot_figure_7(inputs, outputs, threads, config, params, wildcards):
    #load data
    C_cal_par = inputs["C_cal_par"]
    C_cal_root, ext = path.splitext(C_cal_par)
    C_cal = mat.coherencyMatrix(C_cal_root, C_cal_par, gamma=True, bistatic=True,
                                basis='lexicographic').lexicographic_to_pauli()
    #load map
    map_ds = gdal.Open(inputs.map)
    #Create spatial reference
    ref = osr.SpatialReference()
    ref.ImportFromWkt(map_ds.GetProjection())

    #create figure
    plt.style.use(inputs['style'])
    fig_w, fig_h = plt.rcParams['figure.figsize']
    f, ax = plt.subplots(subplot_kw={'projection': ccrs.epsg(ref.GetAuthorityCode('PROJCS'))}, figsize=(2 * fig_w,  2 * fig_w))

    #Geocode
    LUT = geo.GeocodingTable(inputs.dem_seg_par, inputs.LUT)
    C_cal_geo = LUT.geocode_data(C_cal)
    C_cal_geo = C_cal_geo.boxcar_filter([3,3])
    C_cal_rgb = C_cal_geo.pauli_image(k=0.2, sf=0.2, common=False, peak=False)
    C_cal_alpha = np.ones(C_cal_rgb[:,:,0].shape)

    C_cal_alpha[np.isnan(C_cal_geo[:,:,0,0])] = 0
    C_cal_rgb = np.dstack((C_cal_rgb,C_cal_alpha))
    #Plot map
    ext1=geo.get_ds_extent(map_ds)
    ext2=LUT.get_extent()
    print(ext1,ext2)
    ax.imshow(map_ds.ReadAsArray().transpose((1, 2, 0)), extent=ext1)
    ax.imshow(C_cal_rgb.transpose(1,0,2), extent=ext2)
    plt.show()

    # # C_cal_geo, x_vec, y_vec, lut, direct_lut = geo.geocode_image(C_cal, 2)
    # #Decimate reflector location
    # ref_dec = np.array([[ref['ridx'], ref['azidx']//C_cal.azimuth_looks] for ref in params['ref']], dtype=int)
    # #Compute locations in radar coordinates
    # ref_gc = vf.bilinear_interpolate(direct_lut.T, ref_dec[:,0],ref_dec[:,1])
    #
    # ax.imshow(C_cal_rgb.transpose([1,0,2])[::,::-1], aspect=1)
    # marker_color = ['orange' if i == config['calibration']['reflector_index'] else 'cyan' for i in range(len(ref_dec))]
    # ax.scatter(C_cal_geo.shape[0] - np.imag(ref_gc), np.real(ref_gc), edgecolors=marker_color, s=25, facecolors='none', linewidths=0.5)
    # f.subplots_adjust(left=0,right=1)
    ax.axis('off')
    f.savefig(outputs[0])


plot_figure_7(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)
