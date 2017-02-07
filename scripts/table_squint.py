import csv

import numpy as np

import pyrat.fileutils.gpri_files as gpf
from pyrat.gpri_utils import calibration as cal


# Return the decimated azimuth position
def az_idx(ds, idx):
    return idx / ds.azimuth_looks


def center_width_from_slice(sl):
    center = (sl.stop + sl.start) // 2
    width = (sl.stop - sl.start)
    return center, width


width = (10, 100)
z = 1000


def compute_table_4(inputs, outputs, threads, config, params, wildcards):
    # Slc parameters
    slc = gpf.par_to_dict(inputs.slc_par)
    # header for csv
    header = ['name', 'a_HH', 'a_VV', 'a_res_HH', 'a_res_VV']
    res = np.zeros((len(params.ref), 2 * 2))
    for idx_chan, chan in enumerate(('HH', 'VV')):
        for idx_proc, proc in enumerate(('', '_desq')):
            raw = inputs[chan + proc]
            raw_par = raw + '_par'
            raw_data = gpf.rawData(raw_par, raw)
            idx_res = np.ravel_multi_index((idx_proc, idx_chan), (2, 2))
            for idx_ref, ref in enumerate(params.ref):
                squint_idx, squint, sq_par, raw_filt = cal.fit_squint(raw_data, slc, params.ref['azidx'],
                                                                  params.ref['ridx'], win=width, z=z)
                # #Slice
                # az_slice = raw_data.azimuth_slice_from_slc_idx(ref['azidx'], width[1])
                # #Construct
                # raw_sl = raw_data[:, az_slice] * 1
                # #Range filter
                # raw_filt = _sig.hilbert(raw_sl.filter_range_spectrum(slc, ref['ridx'], rwin, k=2), axis=0)
                # az_vec = np.arange(-raw_filt.shape[1]//2, raw_filt.shape[1]//2) * raw_data.azspacing
                # squint_idx, win_slice = squint_vec(raw_filt)
                # squint = az_vec[squint_idx[::-1]]
                # #fit squint
                # sq_par = np.polyfit(raw_data.freqvec[win_slice], squint,1)
                res[idx_ref, idx_res] = sq_par[0]
    names = [ref['name'] for ref in params.ref]
    with open(outputs[0], 'w+') as of:
        writer = csv.writer(of, delimiter=',')
        writer.writerow(header)
        for name, row in zip(names, res):
            out_row = list(row)
            out_row.insert(0, name)
            writer.writerow(out_row)
        mean = list(np.mean(res, axis=0))
        mean.insert(0, 'mean')
        writer.writerow(mean)
    print(res)


compute_table_4(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
                snakemake.wildcards)
