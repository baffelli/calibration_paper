import os.path as path

import matplotlib.pyplot as plt
import numpy as np
import pyrat.core.corefun as cf
import pyrat.core.matrices as mat

import pyrat.visualization.visfun as vf

import pyrat.geo.geofun as geo



# Return the decimated azimuth position
def az_idx(ds, idx):
    print(ds.azimuth_looks)
    return idx / ds.azimuth_looks


def plot_figure_7(inputs, outputs, threads, config, params, wildcards):
    plt.style.use(inputs['style'])
    fig_w, fig_h = plt.rcParams['figure.figsize']
    f, ax = plt.subplots(1, 1, figsize=(fig_w, fig_w ))
    C_cal_par = inputs["C_cal_par"]
    C_cal_root, ext = path.splitext(C_cal_par)
    C_cal = mat.coherencyMatrix(C_cal_root, C_cal_par, gamma=True, bistatic=True,
                                basis='lexicographic').lexicographic_to_pauli()

    C_cal_geo, x_vec, y_vec, lut = geo.geocode_image(C_cal, 2)
    C_cal_rgb = C_cal_geo.pauli_image(k=0.2, sf=0.3, peak=False, perc=[5, 99.9])

    ax.imshow(C_cal_rgb)
    ax.axis('off')
    f.savefig(outputs[0])


plot_figure_7(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)
