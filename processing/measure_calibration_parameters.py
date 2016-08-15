#!/usr/bin/python

import sys, os
import numpy as _np
import argparse
import scipy as _sp
import scipy.signal as _sig
import pyrat.fileutils.gpri_files as _gpf
import pyrat.gpri_utils.calibration as _cal
from collections import namedtuple as _nt
import scipy.signal as _sig
import scipy.ndimage as _ndim
import matplotlib.pyplot as plt
import pyrat.visualization.visfun as vf
from collections import OrderedDict as _od
import pyrat.core.matrices as _mat
import pyrat.core.polfun as _pf
import json
import re
import pyrat.core.corefun as cf
ptarg_re =  re.compile(u"ptarg(\d)(\d)")


def calpal(inputs, outputs, threads, config, params, wildcards):
    C_matrix_flat = _mat.coherencyMatrix(params['C_root'], params['C_par'], basis='lexicographic', gamma=True, bistatic=True)
    #The indices are in the raw coordinates, needs to decimate
    ref_pos, ref_coord = config['list_of_reflectors'][config['calibration']['reflector_index']]
    print(ref_pos)
    azidx = int(ref_pos[1]) / int(config['range_compression']['dec'])
    ridx = int(ref_pos[0])
    print(C_matrix_flat.shape)
    C_ptarg, rplot, azplot, mx_pos, ptarg_info = cf.ptarg(C_matrix_flat, ridx, azidx, azwin=5, rwin=4)
    # plist = _gpf.load_plist(inputs['plist'])
    # # print(plist['ridx'])
    # #Load ptargs
    # for pt_path in inputs['ptarg']:
    #     matches = re.search(ptarg_re, pt_path)
    #     i,j = matches.groups()
    #     C_ptarg[:,:,i,j] = _gpf.load_binary(pt_path, 1024)[:,1024:]
    # max_r, max_az = _np.unravel_index(_np.argmax(C_ptarg[:,:,0,0]), C_ptarg.shape[0:2])

    # #Extract the phase value at the pt
    # HHVV_pt = C_matrix_flat[plist['ridx'],plist['azidx'],3,0]
    # imbalance_pt = (_np.abs(C_matrix_flat[plist['ridx'],plist['azidx'],3,3]) / \
    #                _np.abs(C_matrix_flat[plist['ridx'],plist['azidx'],0,0]))**1/4.0
    # plt.figure()
    # plt.hist(_np.angle(HHVV_pt))
    # rgb, a, cm = vf.dismph(C_matrix_flat[:,:,3,0])
    # plt.figure()
    # plt.imshow(_np.abs(C_matrix_flat[:,:,3,0])**0.1,origin='left')
    # plt.scatter(plist['azidx'], plist['ridx'],
    #             c=_np.angle(C_matrix_flat[plist['ridx'],plist['azidx'],0,3]), s=60)
    # plt.show()
    # HH_VV_diff = _np.mean(_np.angle(HHVV_pt))
    # phi_t_1 = HH_VV_diff + phi_r
    #Compute the HVVH imbalance
    av_win = [5,5]#averaging window
    C_matrix_flat_av = C_matrix_flat.boxcar_filter(av_win)
    phi_t, phi_r, f, g, A = _cal.measure_imbalance(C_ptarg[mx_pos[0], mx_pos[1]], C_matrix_flat_av, config['calibration']['reflector_rcs'])
    cal_dict = _od()
    cal_dict['phi_t'] = phi_t
    cal_dict['phi_r'] = phi_r
    cal_dict['f'] = f.real
    cal_dict['g'] = g.real
    cal_dict['A'] = A



    _gpf.dict_to_par(cal_dict, outputs['cal_par'])



calpal(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params, snakemake.wildcards)






