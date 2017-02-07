import string

import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
from matplotlib import gridspec

import pyrat.fileutils.gpri_files as gpf
from pyrat.gpri_utils import calibration as cal


def az_idx(ds, idx):
    return idx / ds.azimuth_looks


def center_width_from_slice(sl):
    center = (sl.stop + sl.start) // 2
    width = (sl.stop - sl.start)
    return center, width


width = (6, 150)
shp = (2, 2)
# range window for filter in pizzels
# Pixels to discard at beginning and end of chirp
z = 2000
xlabel = r' $\theta_{sq}$ (Azimuth rel. to pointing at $f_c$) $[^\circ]$'
ylabel = r'Chirp frequency $f[MHz]$'


def format_freq(raw_par, x, pos):
    rel_freq = (x - raw_par.RF_center_freq) / 1e6  # relative frequency
    sign = '-' if rel_freq < 0 else '+'
    lab = r'$f_c$ {sign:^10} {rel_freq:^1.0f}'.format(sign=sign, rel_freq=np.abs(rel_freq))
    return lab


def plot_figure_12(inputs, outputs, threads, config, params, wildcards):
    # Create figure
    plt.style.use(inputs['style'])
    fig_w, fig_h = plt.rcParams['figure.figsize']
    f = plt.figure(figsize=(fig_w * 2, fig_h * 2))
    # Grid of plots for figures
    gs = gridspec.GridSpec(*shp, height_ratios=[1, 1])
    gs.update(hspace=0.3, wspace=0.4)
    # data to plot
    for idx_chan, chan in enumerate(('HH', 'VV')):
        for idx_proc, proc in enumerate(('', '_desq')):
            raw = inputs[chan + proc]
            raw_par = raw + '_par'
            raw_data = gpf.rawData(raw_par, raw)
            # Read range and azimuth indices
            ridx = params.ref['azidx']
            azidx = params.ref['azidx']
            # Slc parameters
            slc = gpf.par_to_dict(inputs.slc_par)
            # Fit squint
            squint_idx, squint, sq_par, raw_filt = cal.fit_squint(raw_data, slc, params.ref['azidx'], params.ref['ridx'],
                                                              win=width, z=z)
            ax = f.add_subplot(gs[idx_proc, idx_chan])
            az_vec = np.arange(-raw_filt.shape[1] // 2, raw_filt.shape[1] // 2) * raw_data.azspacing
            squint_fit = raw_data.freqvec * sq_par[0] + sq_par[1]
            ext_vec = [az_vec[0], az_vec[-1], raw_data.freqvec[0], raw_data.freqvec[-1], ]
            ax.imshow(np.abs(raw_filt), aspect=1e-8, extent=ext_vec)
            ax.plot(squint_fit, raw_data.freqvec)
            # Text with fit params
            bbox_props = dict(boxstyle="square", fc="white", ec="w", lw=2)
            t = ax.text(0.1, 0.1, "$a=${:.2f} $\\frac{{\circ}}{{GHz}}$".format(sq_par[0] / 1e-9), size=8,
                        bbox=bbox_props, horizontalalignment='left',
                        transform=ax.transAxes)  # set label
            if idx_chan == idx_proc == 0:
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
            fm_func = lambda x, pos: format_freq(raw_data, x, pos)
            fmt = tick.FuncFormatter(fm_func)
            ax.yaxis.set_major_formatter(fmt)
            # Plot index
            current_idx = np.ravel_multi_index((idx_chan, idx_proc), shp)
            label_name = string.ascii_lowercase[current_idx]
            # current_ax = axarr[idx_chan, idx_proc]
            ax.title.set_text(r"({label_name})".format(label_name=label_name))
            f.suptitle(params.ref['name'])
            f.subplots_adjust(hspace=0.15, wspace=0.05)
            # ax.xaxis.set_major_locator(tick.MultipleLocator(0.5))
            plt.show()
            # New figure to plot subplots
            f1, ax1 = plt.subplots(figsize=(fig_w, fig_h))
            ax1.imshow(np.abs(raw_filt).T, aspect=1e7, extent=ext_vec[::-1])
            ax1.plot(raw_data.freqvec, squint_fit)
            ax1.set_xlabel(ylabel)
            ax1.set_ylabel(r'$\theta_{sq}$ [$^\circ$]')
            f1.savefig(outputs['pres_fig'][current_idx])
        f.savefig(outputs['paper_fig'])


plot_figure_12(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
               snakemake.wildcards)
