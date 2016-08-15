import pyrat.core.matrices as mat
import pyrat.core.polfun as pf
import matplotlib.pyplot as plt
import pyrat.fileutils.gpri_files as gpf
import pyrat.visualization.visfun as vf
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
def display_raw(inputs, outputs, threads, config, params, wildcards):
    with plt.style.context(config['style']):
        raw, raw_dict = gpf.load_raw(inputs['raw'] + '_par', inputs['raw'], nchan=1)
        f, ax = plt.subplots()
        sl = slice(int(wildcards['azidx']) - 200, int(wildcards['azidx']) + 200)
        plt.imshow(raw[::,sl].T, cmap='RdBu_r', aspect=20)
        f.savefig(outputs['image'])
display_raw(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params, snakemake.wildcards)