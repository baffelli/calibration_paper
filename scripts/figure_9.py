import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyrat.fileutils.gpri_files as gpf

import pyrat.visualization.visfun as vf
# Return the decimated azimuth position
def az_idx(ds, idx):
    return idx / ds.azimuth_looks


def plot_figure_8(inputs, outputs, threads, config, params, wildcards):
    C_par = gpf.par_to_dict(inputs.C_cal_par)
    sl = [slice(0,1500), Ellipsis, ]
    HHVV = gpf.load_binary(inputs.HHVV_phase, C_par['range_samples'], dtype=gpf.type_mapping['FCOMPLEX'])
    #Create RGB
    mph_dict = {'k': 0.1, 'sf': 1e-2, 'coherence': False, 'peak': True}
    mph, rgb, norm = vf.dismph(HHVV[sl], **mph_dict)  # create rgb image
    pal, ext = vf.dismph_palette(HHVV[sl], **mph_dict)
    plt.style.use(inputs.style)
    fig_w, fig_h = plt.rcParams['figure.figsize']
    f, ax = plt.subplots(1, 1, figsize=(fig_w,fig_h))
    ax.imshow(mph)
    ax.xaxis.set_label_text(r'range samples')
    ax.yaxis.set_label_text(r'azimuth samples')
    plt.show()
    f.savefig(outputs[0])
    # inc = gpf.load_binary(inputs.inc, C_par['range_samples'], dtype=gpf.type_mapping['FLOAT'])
    # residuals = pd.read_csv(inputs.res)
    # phase_rms = np.std(residuals['HH-VV phase imbalance'])
    # phase_mean = np.mean(residuals['HH-VV phase imbalance'])
    # amp_rms = np.std(residuals['HH-VV amplitude imbalance'])
    # amp_mean = np.mean(residuals['HH-VV amplitude imbalance'])
    # inc_ref = np.rad2deg(inc[residuals['range_index'], residuals['azimuth_index'] / C_par['azimuth_looks']])
    # # perform fit
    # phase_coeff, cov = np.polyfit(inc_ref, residuals['HH-VV phase imbalance'], 3, cov=True)  # thir order for phase
    # plt.style.use(inputs.style)
    # fig_w, fig_h = plt.rcParams['figure.figsize']
    # f, (ph_ax, amp_ax) = plt.subplots(2, 1, figsize=(fig_w * 2, 2* fig_h))
    # ph_ax.scatter(inc_ref, residuals['HH-VV phase imbalance'])
    # ph_ax.xaxis.set_label_text(r'incidence angle $\alpha$ [deg]')
    # ph_ax.yaxis.set_label_text(r'residual $\phi_r + \phi_t$ [deg]')
    # bbox_props = dict(boxstyle="square", fc="white", ec="w", lw=2)
    # ph_ax.text(0.2, 0.2, r"RMS {:.2f}$^\circ$, Mean {:.2f}$^\circ$".format(phase_rms, phase_mean),
    #            transform=ph_ax.transAxes, bbox=bbox_props)
    # ph_ax.set_ylim([-20, 20])
    # amp_ax.scatter(inc_ref, residuals['HH-VV amplitude imbalance'])
    # amp_ax.xaxis.set_label_text(r'incidence angle $\alpha$ [deg] ')
    # amp_ax.yaxis.set_label_text(r'residual $f$')
    # amp_ax.text(0.2, 0.2, r"RMS {:.2f}, Mean {:.2f}".format(amp_rms, amp_mean),
    #            transform=amp_ax.transAxes, bbox=bbox_props)
    # f.subplots_adjust(hspace=0.2)


plot_figure_8(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)
