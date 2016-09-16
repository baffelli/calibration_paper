import matplotlib.pyplot as plt
import pyrat.fileutils.gpri_files as gpf
import pyrat.visualization.visfun as vf
import numpy as np
import scipy.signal as sig
import scipy.optimize as opt
import pyrat.gpri_utils.calibration as cal
import re
import matplotlib as mpl
import mpl_toolkits.axes_grid.inset_locator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.patches import ConnectionPatch

def norm_pat(pat):
    return 10 * np.log10(pat / pat.max())
p = re.compile(r"""[H,V]{2}_(?P<freq>[0-9]{3})GHz""")
def measure_pattern(inputs, outputs, threads, config, params, wildcards):
    raw_data, raw_par = gpf.load_raw(inputs['raw_par'], inputs['raw'], nchan=1)
    raw_data = raw_data / 32768.0
    raw_parameters = gpf.rawParameters(inputs['raw_par'], inputs['raw'])
    #Compute the slc parameters using the config information
    args = type('args', (object,), {})()
    args.zero = int(config['range_compression']['z'])
    args.rmin = float(config['range_compression']['rmin'])
    args.rmax = float(config['range_compression']['rmax'])
    args.kbeta = float(config['range_compression']['k'])
    args.d = int(config['range_compression']['dec'])
    raw_parameters.compute_slc_parameters(args)
    start_az_idx = int(wildcards['azidx'])  + raw_parameters.nl_acc
    sl = slice(start_az_idx - params['azwin']/2, start_az_idx  + params['azwin']/2)
    #Relative azimuth vector
    az_vec = np.linspace(-params['azwin']/2 *  raw_parameters.ang_per_tcycle,
                       params['azwin']/2 * raw_parameters.ang_per_tcycle,params['azwin'])
    #Compress data
    slc = np.fft.rfft(raw_data[1:, sl],axis=0)
    fshift = np.zeros_like(slc) + 1
    #Range video phase

    #Position of range index (taking the rmin parameter into account
    ridx = int(wildcards['ridx']) + raw_parameters.ns_min
    #filter only the ranges of interest
    filt = np.zeros(slc.shape[0])
    filt[ridx-params['rwin']/2:ridx+params['rwin']/2] = np.hamming(params['rwin'])
    #Filter the samples around the reflector
    raw_filt = np.fft.irfft(slc  * filt[:, None],axis=0)
    #Compute the envelope
    raw_envelope = sig.hilbert(raw_filt[:, :],axis=0)
    raw_envelope = raw_envelope * raw_envelope[:, az_vec.shape[0]/2].conj()[:,None]
    #Load measured patterns
    pat_name = 'HH' if wildcards['chan'][0:3] == 'AAA' else 'VV'
    measured_patterns = {}
    for pat_path in inputs[pat_name]:
        ma = re.search(p, pat_path)
        freq = int(ma.group('freq'))
        measured_patterns[freq] = np.genfromtxt(pat_path, delimiter=',')
    #Measure pattern
    min_pat = norm_pat(np.abs(raw_filt[0,:])**2)
    mid_pat = norm_pat(np.abs(raw_filt[12345,:])**2)
    f, response_axis = plt.subplots()
    #dictionary that matches the measured patterns with the frequency vector
    freq_dict = {171:1500, 172:raw_parameters.nsamp/2,173:raw_parameters.nsamp - 1500}
    with plt.style.context(config['style']):
        #Main axis with the response
        rgb, c, d = vf.dismph(raw_envelope)
        gs = gridspec.GridSpec(20,2)
        response_axis = plt.subplot(gs[:,0])
        response_axis.imshow(rgb, cmap=c, extent=[az_vec.min(), az_vec.max(),
                                                 raw_parameters.freq_vec.min(), raw_parameters.freq_vec.max(),],
                   aspect=1e-7)
        #create axis divider
        # freq_dict = np.linspace(1200,24999,20)
        # for idx, freq in enumerate(freq_dict):
        #     current_envelope = raw_envelope[freq_dict[idx]]
        #     estimated_pat = norm_pat(np.abs(current_envelope)**2)
        #     # measured_pattern = measured_patterns[freq]
        #     # plot_axis = plt.subplot(gs[idx % 1 + idx, 1])
        #     phase_axis = plt.subplot(gs[idx, 1])
        #     # est_pat_plot  = plot_axis.plot(az_vec, estimated_pat ,'--')
        #     # plot_axis.plot(measured_pattern[:,0],measured_pattern[:,1], color=est_pat_plot[0].get_color())
        #     norm_angle = np.angle(current_envelope * current_envelope[current_envelope.shape[0]/2].conj())
        #     phase_plot = phase_axis.plot(az_vec,np.rad2deg((norm_angle)))
        #     phase_axis.set_ylim(-180,180)
        #     #Plot line in response axis
        #     current_frequency = raw_parameters.freq_vec[freq_dict[idx]]
        #     response_axis.axhline(y=current_frequency,color=phase_plot[0].get_color())
        #     #plot connector
        #     # freq_az = (current_frequency, 0)
        #     # az_amp = (0,0)
        #     # con = ConnectionPatch(xyA=freq_az, xyB=az_amp, coordsA='data',
        #     #                       coordsB='data', axesA=response_axis, axesB=plot_axis)
        #     # plot_axis.add_artist(con)
    plt.show()


measure_pattern(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params, snakemake.wildcards)

