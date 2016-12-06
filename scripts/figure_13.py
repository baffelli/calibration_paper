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
import scipy.interpolate as interp




width = (10, 150)
shp = (2,2)
#range window for filter in pizzels
rwin = 4

xlabel = r' $\theta_{sq}$ (Azimuth rel. to pointing at $f_c$) $[\circ]$'
ylabel = r'Chirp frequency $f[MHz]$'

interp_funct = lambda x,y: interp.interp1d(x,y/2, kind='slinear')

def plot_figure_13(inputs, outputs, threads, config, params, wildcards):
    # Create figure
    plt.style.use(inputs['style'])
    fig_w, fig_h = plt.rcParams['figure.figsize']
    f = plt.figure(figsize=(2 * fig_w, 2 * fig_h))
    ax = f.add_subplot(111)
    #Load csv
    H_pat = np.genfromtxt(inputs['H_pat'], delimiter=',')
    V_pat = np.genfromtxt(inputs['V_pat'], delimiter=',')
    #Load coregistration parameters
    coreg_par = gpf.par_to_dict(inputs['coreg_par'])
    slc_par = gpf.par_to_dict(inputs['slc_par'] + '.par')
    #Compute shift in samples
    shift = coreg_par.azimuth_offset_polynomial[0] * slc_par.GPRI_az_angle_step
    print(shift)
    #Construct interpolants
    V_interp = interp_funct(V_pat[:,0], V_pat[:,1])
    H_interp = interp_funct(H_pat[:, 0], H_pat[:, 1])
    az_vec = np.linspace(-0.5,0.5)
    # V_pat_shift = np.interp(V_pat[:,0], V_pat[:, 0] + shift,V_pat[:,1], )
    ax.plot(az_vec,H_interp(az_vec), label='H antenna')
    ax.plot(az_vec, V_interp(az_vec), label='V antenna')
    ax.plot(az_vec + shift, V_interp(az_vec), label='V antenna, shifted')
    gain = H_interp(0) - V_interp(shift)
    print(gain)
    # ax.plot(V_pat[:, 0], V_pat_shift, label='V antenna')
    ax.xaxis.set_label_text(r'Angle from maximum [$\circ$]')
    ax.yaxis.set_label_text(r'Relative One-way Pattern [dB]')
    ax.annotate('', xy=(0, H_interp(0)), xytext=(0, V_interp(0) - gain),
                xycoords='data', textcoords='data', arrowprops=dict(arrowstyle="|-|",
                                                                    ec="k",
                                                                    lw=0.4,
                                                                    shrinkA = 0,
                                                                    shrinkB=0,
                                                                    ))
    ax.text(0 - 0.1, (H_interp(0) - gain)/2,
            "loss: {HV:.2f} dB".format(HV=gain), bbox=dict(boxstyle="square", fc="w", ec="k", lw=0.4))
    plt.legend()
    f.savefig(outputs[0])



plot_figure_13(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)
