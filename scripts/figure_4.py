import matplotlib.pyplot as plt
import numpy as np
import pyrat.fileutils.gpri_files as gpf
import pyrat.core.corefun as cf
import pyrat.core.polfun as  pf

def plot_figure_4(inputs, outputs, threads, config, params, wildcards):
    plt.style.use(inputs['style'])
    HV = gpf.gammaDataset(inputs.C_HV_old + '.par', inputs.C_HV_old)
    HV_corr = gpf.gammaDataset(inputs.C_HV_new + '.par', inputs.C_HV_new)
    ridx = params.ridx
    azidx = params.azidx
    sw = (2, 30)
    azwin = 16
    ptarg_zoom, rplot, azplot, mx_idx_zoom, res_dict, r_vec, az_vec = cf.ptarg(HV, ridx, azidx, sw=sw, azwin=azwin)
    ptarg_zoom_corr, rplot_corr, azplot_corr, mx_idx_zoom_corr, res_dict, r_vec_corr, az_vec_corr = cf.ptarg(HV_corr, ridx, azidx, sw=sw,
                                                                                    azwin=azwin)
    HV_gain = cf.dB(ptarg_zoom_corr[mx_idx_zoom_corr] / ptarg_zoom[mx_idx_zoom])
    x_shift = 50
    fig_w, fig_h = plt.rcParams['figure.figsize']
    f, ax = plt.subplots(1, 1, figsize=(2 * fig_w, 2 * fig_h))
    bbox_props = dict(boxstyle="square", fc="white", ec="w", lw=2)
    norm_plot = plt.plot(az_vec, cf.dB(np.abs(azplot)), label='No shift')
    corr_plot = plt.plot(az_vec_corr, cf.dB(np.abs(azplot_corr)), label='Optimal shift')
    plt.axhline(cf.dB(np.abs(ptarg_zoom[mx_idx_zoom])), color=norm_plot[0].get_color(), linestyle='--')
    plt.axhline(cf.dB(np.abs(ptarg_zoom_corr[mx_idx_zoom_corr])), color=corr_plot[0].get_color(), linestyle='--')
    ax.arrow(az_vec[mx_idx_zoom[1] + x_shift], cf.dB(np.abs(ptarg_zoom[mx_idx_zoom])), 0, HV_gain)
    ax.text(az_vec[mx_idx_zoom[1] + x_shift + 18], cf.dB(np.abs(ptarg_zoom[mx_idx_zoom])) + HV_gain / 2,
             "gain: {HV:.2f} dB".format(HV=HV_gain), bbox=bbox_props)
    ax.xaxis.set_label_text(r'azimuth samples')
    ax.yaxis.set_label_text(r'HV power [dB]')
    ax.set_ylim([0, 35])
    ax.grid(True)
    plt.legend(loc='lower left',frameon=True)
    f.subplots_adjust(bottom=0.15)
    f.savefig(outputs[0], pad_inches=0.1)

plot_figure_4(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)