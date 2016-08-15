import pyrat.core.matrices as mat
import pyrat.core.polfun as pf
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pyrat.fileutils.gpri_files as gpf
import pyrat.visualization.visfun as vf
import numpy as np
def plot_phase(inputs, outputs, threads, config, params, wildcards):
    dem_seg_par = gpf.par_to_dict(inputs['dem_seg_par'])
    ifgram = gpf.load_binary(inputs['data'], dem_seg_par['width'], dtype=gpf.type_mapping['FCOMPLEX']).T
    print(ifgram.shape)
    #Create rgb image
    mph, rgb, norm = vf.dismph(ifgram, k=0.2, sf=1)
    mph[mph==0] = 1
    with plt.style.context(inputs['style']):
        #This one has to be square
        w,h = plt.rcParams['figure.figsize']
        f, ax = plt.subplots(figsize=(w,w))
        plt.imshow(mph, interpolation='none')
        plt.tight_layout()
        ax.axis('off')
        f.savefig(outputs['image'])


plot_phase(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params, snakemake.wildcards)