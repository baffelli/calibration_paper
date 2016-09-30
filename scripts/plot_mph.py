import pyrat.core.matrices as mat
import pyrat.core.polfun as pf
import matplotlib.pyplot as plt
import pyrat.core.corefun as cf
import pyrat.fileutils.gpri_files as gpf
import pyrat.visualization.visfun as vf
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import scipy.ndimage as ndim
def plot_mph(inputs, outputs, threads, config, params, wildcards):
    with plt.style.context(config['style']):
        slc = gpf.gammaDataset(inputs['slc'] + '.par', inputs['slc'] )
        r_idx = int(wildcards['ridx'])
        az_idx = int(wildcards['azidx'])
        #Compute the point target response
        osf = 64
        ptarg_zoom, r_plot, az_plot, mx_pos, res_dict = cf.ptarg(slc,r_idx,az_idx,azwin=100, rwin=20)
        az_vec = np.arange(-ptarg_zoom.shape[1]/2,ptarg_zoom.shape[1]/2) * slc.GPRI_az_angle_step[0] / osf
        r_vec = np.arange(-ptarg_zoom.shape[0] / 2, ptarg_zoom.shape[0] / 2) * slc.range_pixel_spacing[0] / osf
        mph,rgb, norm = vf.dismph(ptarg_zoom.T , k=0.4, sf=0.1, coherence=False)
        w,h = plt.rcParams['figure.figsize']
        #Plot range and azimuth responses too
        f, (mph_ax) =  plt.subplots(1,1, figsize=(w,w))
        mappable = mph_ax.imshow(mph[:, :], cmap=rgb, extent=[r_vec.min(), r_vec.max(), az_vec.min(),az_vec.max()], aspect=10.0, interpolation='none')
        # #Create axes for range and azimuth respoinse
        # divider = make_axes_locatable(mph_ax)
        plt.figtext(0.15,0.95, r"\begin{{tabular}}{{l}}Range resolution: {r_res:.3f} m\\Azimuth resolution:"
                             " {az_res:.3f}$^\circ$\end{{tabular}}".format(r_res=res_dict['range_resolution'][0],
                                                                            az_res=res_dict['azimuth_resolution'][0]), fontsize=5)
        mph_ax.set_xlabel(r'range samples [m]')
        mph_ax.set_ylabel(r'azimuth samples $\theta$ [deg]')
        ticks = [0, 0.25, 0.5, 0.75, 1]
        labels = ['-180', '-90', '0', '90', '180']
        #Make space for colorbar
        cbar = f.colorbar(mappable, label=r'Phase [deg]', orientation='vertical', ticks=ticks, shrink=0.5, pad=0.1)
        cbar.set_ticks(ticks)
        cbar.set_ticklabels(labels)
        cbar.update_ticks()
        plt.grid()
        f.savefig(outputs['plot'])

plot_mph(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params, snakemake.wildcards)