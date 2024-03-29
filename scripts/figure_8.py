import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyrat.fileutils.gpri_files as gpf

import matplotlib.ticker as tick

# Return the decimated azimuth position
def az_idx(ds, idx):
    return idx / ds.azimuth_looks

x_box = 0.1
y_box = 0.2

def plot_figure_8(inputs, outputs, threads, config, params, wildcards):
    C_par = gpf.par_to_dict(inputs.C_cal_par)
    inc = gpf.load_binary(inputs.inc, C_par['range_samples'], dtype=gpf.type_mapping['FLOAT'])
    residuals = pd.read_csv(inputs.res)
    phase_rms = np.std(residuals['HH-VV phase imbalance'])
    phase_mean = np.mean(residuals['HH-VV phase imbalance'])
    amp_rms = np.std(residuals['HH-VV amplitude imbalance'])
    amp_mean = np.mean(residuals['HH-VV amplitude imbalance'])
    inc_ref = np.rad2deg(inc[residuals['range_index'], residuals['azimuth_index'] // C_par['GPRI_decimation_factor']])
    # perform fit
    phase_coeff, cov = np.polyfit(inc_ref, residuals['HH-VV phase imbalance'], 3, cov=True)  # thir order for phase
    plt.style.use(inputs.style)
    fig_w, fig_h = plt.rcParams['figure.figsize']
    xlabel = r'Local incidence angle $[^\circ]$'
    f, (ph_ax, amp_ax) = plt.subplots(2, 1, figsize=(fig_w, fig_h))
    ph_ax.scatter(inc_ref, residuals['HH-VV phase imbalance'])
    # ph_ax.xaxis.set_label_text(xlabel)
    ph_ax.yaxis.set_label_text(r'Calibrated $\phi_r + \phi_t [^\circ]$')
    bbox_props = dict(boxstyle="square", fc="white", ec="black", lw=0.5)
    ph_ax.text(x_box, y_box, r"RMS {:.2f}$^\circ$, Mean {:.2f}$^\circ$".format(phase_rms, phase_mean),
               transform=ph_ax.transAxes, bbox=bbox_props)
    ph_ax.set_ylim([-20, 20])
    ph_ax.yaxis.set_major_locator(tick.MultipleLocator(10))
    amp_ax.scatter(inc_ref, residuals['HH-VV amplitude imbalance'])
    amp_ax.xaxis.set_label_text(xlabel)
    amp_ax.yaxis.set_label_text(r'Calibrated $f$')
    amp_ax.text(x_box, y_box, r"RMS {:.2f}, Mean {:.2f}".format(amp_rms, amp_mean),
               transform=amp_ax.transAxes, bbox=bbox_props)
    amp_ax.set_ylim([0.9,1.2])
    f.subplots_adjust(hspace=0.35, right=0.85)
    f.savefig(outputs[0])


plot_figure_8(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)
