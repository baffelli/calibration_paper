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

    #load dem par
    dem_par = gpf.par_to_dict(inputs['dem_seg_par'])
    cov_par = gpf.par_to_dict(inputs['cov_par'])
    #load list of reflectors
    with open(inputs['reflectors']) as refl_list:
        reflectors = json.load(refl_list)['list_of_reflectors']
    #Decimate the azimuth indices
    ridx = [int(ref[0][0]) for ref in reflectors]
    azidx = [int(ref[0][1])/int(config['range_compression']['dec']) for ref in reflectors]
    #load coherency matrix
    C_mat = np.array([gpf.load_binary(params['C_root'] + ".c{i}{j}_gct".format(i=i,j=j),
                                      dem_par['width']) for i in range(4) for j in range(4)]).transpose([2,1,0])
    C_mat = mat.coherencyMatrix(C_mat.reshape((C_mat.shape[0],C_mat.shape[1],4,4)),basis='lexicographic', bistatic=True).to_monostatic().lexicographic_to_pauli()
    rgb = vf.mask_zeros(C_mat.pauli_image(k=0.1, sf=2))
    # DS = gdal.Open(inputs['swiss_image'])
    # img = DS.ReadAsArray().transpose([1,2,0])
    # si_ext = geo.get_ds_extent(DS)
    # dem_ext = geo.get_dem_extent(dem_par)
    # #load LUT
    # LUT = gpf.load_binary(inputs['lut'], cov_par['range_samples'],dtype=gpf.type_mapping['FCOMPLEX'])
    # conv_coord = []
    # for ridx, azidx in zip(ridx, azidx):
    #     coord = LUT[ridx, azidx]
    #     conv_coord.append((coord.imag , coord.real ))
    # refidx = 2
    # sl = [slice(conv_coord[refidx][idx] - 50,conv_coord[refidx][idx] + 50) for idx in [0,1]]
    # conv_coord = np.array(conv_coord)
    #This one has to be square
    with plt.style.context(config['style']):
        w,h = plt.rcParams['figure.figsize']
        f, ax = plt.subplots(figsize=(w,w))
        ax.imshow(rgb, interpolation='none')
        ax.axis('off')
        f.tight_layout()
        f.savefig(outputs['image'])
        # plt.show()


plot_reflectors(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params, snakemake.wildcards)