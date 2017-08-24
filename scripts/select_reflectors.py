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

import pyrat.geo.utils as ut

import pyrat.geo.transforms as trasf

def select_reflectors(inputs, outputs, threads, config, params, wildcards):
    #load map
    map_ds = gdal.Open(inputs.map)
    dem_seg = gdal.Open(inputs.dem_seg)
    #Load geocoding table
    LUT = geo.GeocodingTable(inputs.dem_seg_par, inputs.LUT, inputs.ref_mli_par, inputs.LUT1)
    ext1=geo.get_ds_extent(map_ds)
    or_pt = [1500,300]
    pt1 =LUT.radar_coord_to_geo_coord(or_pt)

    f, ax = plt.subplots()
    ax.imshow(map_ds.ReadAsArray().transpose((1, 2, 0)), extent=ext1)
    ax.imshow(dem_seg.ReadAsArray().real, extent=geo.get_ds_extent(dem_seg))
    ax.plot(pt1[:,0], pt1[:,1], marker='o')
    plt.show()
    # print(dir(LUT.radar_idx_to_dem_idx_t))
    # radar_idx = LUT.radar_coord_to_geo_coord([10,10])
    # print(radar_idx)
    # #Get geo transform
    # gt = geo.get_geotransform(gpf.par_to_dict(inputs.dem_seg_par))
    # # a = trasf.LutAxes(gt=gt, lut=np.array(LUT))
    # # a.imshow(np.abs(slc))
    # f = plt.axes(projection=trasf.LutProjection(gt=gt, lut=np.array(LUT)))


    # LUT = geo.GeocodingTable(inputs.dem_seg_par, inputs.LUT)
    # plt.subplot(2,1,1)
    # plt.imshow(np.real(LUT.lut.lut))
    # plt.subplot(2, 1, 2)
    # plt.imshow(np.real(LUT.inverse_lut.lut))
    # plt.show()
    # pt = [2500,363]
    # radar_index = LUT.dem_index_to_radar_index(pt)
    # print(radar_index)
    # dem_index = LUT.radar_index_to_dem_index(radar_index[0])
    # print(dem_index,pt)
    # #load map
    # map_ds = gdal.Open(inputs.map)
    # #Create spatial reference
    # ref = osr.SpatialReference()
    # ref.ImportFromWkt(map_ds.GetProjection())
    # ut.GeocodedVisualization(slc, LUT, map_ds)
    # plt.show()




select_reflectors(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)
