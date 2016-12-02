import matplotlib.pyplot as plt
import numpy as np
import pyrat.core.corefun as cf
import pyrat.fileutils.gpri_files as gpf


def plot_figure_4(inputs, outputs, threads, config, params, wildcards):
    plt.style.use(inputs['style'])
    HV = gpf.gammaDataset(inputs.C_HV_old + '.par', inputs.C_HV_old)
    HV_corr = gpf.gammaDataset(inputs.C_HV_new + '.par', inputs.C_HV_new)
    ridx = params.ref['ridx']
    azidx = params.ref['azidx']
    sw = (2, 200)
    azwin = 100
    ptarg_zoom, rplot, azplot, mx_idx_zoom, res_dict, r_vec, az_vec = cf.ptarg(np.abs(HV) ** 2, ridx, azidx, sw=sw,
                                                                               azwin=azwin)
    ptarg_zoom_corr, rplot_corr, azplot_corr, mx_idx_zoom_corr, res_dict, r_vec_corr, az_vec_corr = cf.ptarg(
        np.abs(HV_corr) ** 2, ridx, azidx, sw=sw,
        azwin=azwin)
    HV_gain = cf.dB(ptarg_zoom_corr[mx_idx_zoom_corr] / ptarg_zoom[mx_idx_zoom])
    x_shift = 200
    fig_w, fig_h = plt.rcParams['figure.figsize']
    f, ax = plt.subplots(1, 1, figsize=(2 * fig_w, 2 * fig_h))
    bbox_props = dict(boxstyle="square", fc="w", ec="k", lw=0.4)
    norm_plot = plt.plot(az_vec, cf.dB(np.abs(azplot)), label='No shift')
    corr_plot = plt.plot(az_vec_corr, cf.dB(np.abs(azplot_corr)), label='Optimal shift')
    # lines at the maximum
    max = cf.dB(np.abs(ptarg_zoom[mx_idx_zoom]))
    max_corr = cf.dB(np.abs(ptarg_zoom_corr[mx_idx_zoom_corr]))
    max_az = az_vec[mx_idx_zoom[1] + x_shift]
    plt.axhline(max, color=norm_plot[0].get_color(), linestyle='--')
    plt.axhline(max_corr, color=corr_plot[0].get_color(), linestyle='--')
    # Add arrow with gain
    ax.annotate('', xy=(max_az, max), xytext=(max_az, max_corr),
                xycoords='data', textcoords='data', arrowprops=dict(arrowstyle="|-|",
                                                                    ec="k",
                                                                    lw=0.4,
                                                                    shrinkA = 0,
                                                                    shrinkB=0,
                                                                    ))
    ax.text(az_vec[mx_idx_zoom[1] + x_shift + 80], (max_corr + max)/2,
            "gain: {HV:.2f} dB".format(HV=HV_gain), bbox=bbox_props)
    ax.xaxis.set_label_text(r'azimuth [$\circ$]')
    ax.yaxis.set_label_text(r'HV power [dB]')
    ax.set_ylim([-20, 35])
    ax.grid(True)
    plt.legend(loc='upper left', frameon=True)
    f.subplots_adjust(bottom=0.15)
    f.savefig(outputs[0], pad_inches=0.1)


plot_figure_4(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)
