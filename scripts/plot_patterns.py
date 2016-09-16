import pyrat.core.matrices as mat
import pyrat.core.polfun as pf
import matplotlib.pyplot as plt
import pyrat.fileutils.gpri_files as gpf
import pyrat.visualization.visfun as vf
import numpy as np
import re
p = re.compile(r"""[H,V]{2}_(?P<freq>[0-9]{3})GHz""")
def plot_pat(inputs, outputs, threads, config, params, wildcards):
    with plt.style.context(config['style']):
        f, ax = plt.subplots()
        for HH_path, VV_path in zip(inputs['HH'], inputs['VV']):
            ma = re.search(p, HH_path)
            freq = int(ma.group('freq'))
            HH = np.genfromtxt(HH_path, delimiter=',')
            VV = np.genfromtxt(VV_path, delimiter=',')
            HH_plot = ax.plot(HH[:,0], HH[:,1], label='HH@{freq} GHz'.format(freq = freq/10 ))
            ax.plot(VV[:,0], VV[:,1], color=HH_plot[0].get_color(),ls="--", label='VV@{freq} GHz'.format(freq = freq/10 ))
        ax.set_xlabel(r'azimuth [deg]')
        ax.set_ylabel(r'antenna pattern [dB]')
        ax.legend()
        ax.set_ylim([-40,0])
        ax.set_xlim([-4.5, 4.5])
        f.savefig(outputs['plot'])
plot_pat(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params, snakemake.wildcards)

