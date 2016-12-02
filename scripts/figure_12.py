import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
import pyrat.core.corefun as cf
import pyrat.diff.intfun as _ifun
import pyrat.fileutils.gpri_files as gpf
import pyrat.visualization.visfun as vf
from matplotlib import gridspec
import matplotlib.ticker as tick
import string
import scipy.signal as _sig


# Return the decimated azimuth position
def az_idx(ds, idx):
    return idx / ds.azimuth_looks

def center_width_from_slice(sl):
    center = (sl.stop + sl.start)//2
    width = (sl.stop - sl.start)
    return center, width


def squint_vec(rawdata, z=1000):
    win_slice = slice(z, rawdata.shape[0] - z)
    rawdata_sl = rawdata[win_slice,:]#window the edges
    max_idx = np.argmax(np.abs(rawdata_sl), axis=1)
    return max_idx, win_slice


width = (10, 100)
shp = (2,2)
#range window for filter in pizzels
rwin = 5

def format_freq(raw_par, x, pos):
    rel_freq = (x - raw_par.RF_center_freq) / 1e6#relative frequency
    sign = '-' if rel_freq < 0 else '+'
    lab = r'$f_c$ {sign:^10} {rel_freq:^1.0f}'.format(sign= sign, rel_freq=np.abs(rel_freq))
    return lab

def plot_figure_12(inputs, outputs, threads, config, params, wildcards):
    # Create figure
    plt.style.use(inputs['style'])
    fig_w, fig_h = plt.rcParams['figure.figsize']
    f = plt.figure(figsize=(2 * fig_w, 2 * fig_h))
    # Grid of plots for figures
    gs = gridspec.GridSpec(*shp, height_ratios=[1, 1])
    gs.update(hspace=0.3, wspace=0.4)
    #data to plot
    for idx_chan, chan in enumerate(('HH', 'VV')):
        for idx_proc, proc in enumerate(('', '_desq')):
            raw = inputs[chan + proc]
            raw_par = raw + '_par'
            raw_data = gpf.rawData(raw_par, raw)
            #Slc parameters
            slc = gpf.par_to_dict(inputs.slc_par)
            #Slice
            az_slice = raw_data.azimuth_slice(params.ref['azidx'], width[1])
            #Construct
            raw_sl = raw_data[:, az_slice] * 1
            #Range filter
            raw_filt = _sig.hilbert(raw_sl.filter_range_spectrum(slc, params.ref['ridx'], rwin), axis=0)
            ax = f.add_subplot(gs[idx_proc, idx_chan])
            az_vec = np.arange(-raw_filt.shape[1]//2, raw_filt.shape[1]//2) * raw_data.azspacing
            squint_idx, win_slice = squint_vec(raw_filt)
            squint = az_vec[squint_idx[::-1]]
            #fit squint
            sq_par = np.polyfit(raw_data.freqvec[win_slice], squint,1)
            squint_fit = raw_data.freqvec * sq_par[0] + sq_par[1]
            print(sq_par)
            ax.imshow(np.abs(raw_filt), aspect=1e-8, extent = [az_vec[0], az_vec[-1], raw_data.freqvec[0], raw_data.freqvec[-1],])
            ax.plot(squint_fit, raw_data.freqvec)
            #Text with fit params
            bbox_props = dict(boxstyle="square", fc="white", ec="w", lw=2)
            t = ax.text(0.1, 0.1, "$a=${:.2f} $\\frac{{\circ}}{{GHz}}$".format(sq_par[0]/1e-9), size=8, bbox=bbox_props, horizontalalignment='left',
                                transform=ax.transAxes)  # set label
            if idx_chan == idx_proc == 0:
                ax.set_xlabel(r' $\theta_{sq}$ (Azimuth from pointing at $f_c$) $[\circ]$')
                ax.set_ylabel(r'Chirp frequency $f[MHz]$')
            fm_func = lambda x, pos: format_freq(raw_data, x, pos)
            fmt = tick.FuncFormatter(fm_func)
            ax.yaxis.set_major_formatter(fmt)
            #Plot index
            current_idx = np.ravel_multi_index((idx_chan , idx_proc), shp)
            label_name = string.ascii_lowercase[current_idx]
            # current_ax = axarr[idx_chan, idx_proc]
            ax.title.set_text(r"({label_name})".format(label_name=label_name))
            f.suptitle(params.ref['name'])
            # ax.xaxis.set_major_locator(tick.MultipleLocator(0.5))
            f.savefig(outputs[0])
    plt.show()



plot_figure_12(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)
