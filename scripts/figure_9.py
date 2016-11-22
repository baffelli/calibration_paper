import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import pyrat.core.corefun as cf
import pyrat.fileutils.gpri_files as gpf
import pyrat.visualization.visfun as vf
import pyrat.diff.intfun as _ifun
import numpy as np
from matplotlib import gridspec


# Return the decimated azimuth position
def az_idx(ds, idx):
    return idx / ds.azimuth_looks


def plot_figure_9(inputs, outputs, threads, config, params, wildcards):
    print(inputs.HHVV_phase)
    HHVV = gpf.gammaDataset(inputs.C_cal_par, inputs.HHVV_phase[0], dtype=gpf.type_mapping['FCOMPLEX'])
    HH = gpf.gammaDataset(inputs.C_cal_par, inputs.HHVV_phase[1], dtype=gpf.type_mapping['FCOMPLEX'])
    VV = gpf.gammaDataset(inputs.C_cal_par, inputs.HHVV_phase[2], dtype=gpf.type_mapping['FCOMPLEX'])
    win = [3,2]#multilooking window
    HHVV = _ifun.estimate_coherence(HHVV, HH, VV, win, discard=True)
    mli = cf.smooth(VV, win , discard=True)
    az_vec = HHVV.az_vec
    r_vec = HHVV.r_vec
    # Create RGB
    mph_dict = {'k': 0.07, 'sf': 0.9, 'coherence': True, 'peak': False, 'mli':mli, 'coherence_threshold':0.45, 'coherence_slope':25}
    mph, rgb, norm = vf.dismph(HHVV, **mph_dict)  # create rgb image
    pal, ext = vf.dismph_palette(HHVV, **mph_dict)
    print(ext)

    plt.style.use(inputs['style'])
    fig_w, fig_h = plt.rcParams['figure.figsize']
    f =  plt.figure(figsize=(2 * fig_w, 2 * fig_h))
    gs = gridspec.GridSpec(*(2, 4), height_ratios=[1,0.2])
    gs.update(hspace=0.3, wspace=0.0)
    im_ax = f.add_subplot(gs[0,::])
    aspect = fig_h / fig_w
    im_ax.imshow(mph, extent=[az_vec[0], az_vec[-1], r_vec[-1], r_vec[1]], aspect=1/25, origin='upper')
    im_ax.yaxis.set_label_text(r'range [m]')
    im_ax.xaxis.set_label_text(r'azimuth [$^\circ$]')
    im_ax.yaxis.set_major_locator(tick.MultipleLocator(500))
    im_ax.xaxis.set_major_locator(tick.MultipleLocator(20))
    #plot palette
    pal_ax = f.add_subplot(gs[-1,1])
    pal_ax.imshow(pal, extent=ext, aspect=1e4)
    pal_ax.set_ylabel(r'Phase')
    pal_ax.set_xlabel(r'Intensity')
    pal_ax.grid(b=False)
    pal_ax.set_yticks([-np.pi, 0, np.pi])
    pal_ax.set_yticklabels([r"$-\pi$", r"$0$", r"$\pi$"])
    ext_list = [ext[0], (ext[0] + ext[1])/2, ext[1]]
    pal_ax.set_xticks([])
    # pal_ax.set_xticklabels([r"{{val:.2f}}".format(val=val) for val in ext_list])
    #Plot coherence
    coh_ax = f.add_subplot(gs[-1, 2])
    c = np.linspace(0,1)
    c_scale = vf.scale_coherence(c, threshold=mph_dict['coherence_threshold'], slope=mph_dict['coherence_slope'])
    coh_ax.plot(c,c_scale)
    coh_ax.xaxis.set_label_text(r'Coherence')
    coh_ax.yaxis.set_label_text(r'Saturation')
    coh_ax.set_aspect(1)
    f.savefig(outputs[0])


plot_figure_9(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)
