import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import pyrat.core.corefun as cf
import pyrat.fileutils.gpri_files as gpf
import pyrat.visualization.visfun as vf
import pyrat.diff.intfun as _ifun
import numpy as np

# Return the decimated azimuth position
def az_idx(ds, idx):
    return idx / ds.azimuth_looks


def plot_figure_8(inputs, outputs, threads, config, params, wildcards):
    print(inputs.HHVV_phase)
    HHVV = gpf.gammaDataset(inputs.C_cal_par, inputs.HHVV_phase[0], dtype=gpf.type_mapping['FCOMPLEX'])
    HH = gpf.gammaDataset(inputs.C_cal_par, inputs.HHVV_phase[1], dtype=gpf.type_mapping['FCOMPLEX'])
    VV = gpf.gammaDataset(inputs.C_cal_par, inputs.HHVV_phase[2], dtype=gpf.type_mapping['FCOMPLEX'])
    HHVV = _ifun.estimate_coherence(HHVV, HH, VV, [5,5])
    az_vec = HHVV.az_vec
    r_vec = HHVV.r_vec
    # Create RGB
    mph_dict = {'k': 0.2, 'sf': 1.2, 'coherence': True, 'peak': False, 'mli':HH, 'coherence_threshold':0.6}
    mph, rgb, norm = vf.dismph(HHVV, **mph_dict)  # create rgb image
    # pal, ext = vf.dismph_palette(HHVV, **mph_dict)


    plt.style.use(inputs['style'])
    fig_w, fig_h = plt.rcParams['figure.figsize']
    f, im_ax = plt.subplots(1, 1, figsize=(2 * fig_w, 2 * fig_h), sharey=True)
    im_ax.imshow(mph, extent=[az_vec[0], az_vec[-1], r_vec[-1], r_vec[1]], aspect=1 / 10, origin='upper')
    im_ax.yaxis.set_label_text(r'range [m]')
    im_ax.xaxis.set_label_text(r'azimuth [$^\circ$]')
    im_ax.yaxis.set_major_locator(tick.MultipleLocator(500))
    im_ax.xaxis.set_major_locator(tick.MultipleLocator(20))
    f.savefig(outputs[0])


plot_figure_8(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)
