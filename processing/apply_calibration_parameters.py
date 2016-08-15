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


def apply_cal(inputs, outputs, threads, config, params, wildcards):
    C_matrix_flat = _mat.coherencyMatrix(params['C_input_root'], params['C_input_root'] + '.par', basis='lexicographic', gamma=True, bistatic=True)
    cal_dict = _gpf.par_to_dict(inputs['cal_par'])
    phi_t = cal_dict['phi_t']
    phi_r = cal_dict['phi_r']
    f = cal_dict['f']
    g = cal_dict['g']
    A = cal_dict['A']
    D = _cal.distortion_matrix(phi_t, phi_r, f, g)
    D_inv = _np.diag(1/_np.diag(D))
    C_cal = A**4 * C_matrix_flat.transform(D_inv,D_inv.T.conj())
    C_cal.to_gamma(params['C_output_root'], bistatic=True)




apply_cal(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params, snakemake.wildcards)






