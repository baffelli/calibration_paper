import os.path as path

import matplotlib.pyplot as plt
import numpy as np
import pyrat.core.corefun as cf
import pyrat.core.matrices as mat
import pyrat.visualization.visfun as vf
from cycler import cycler
from mpl_toolkits.axes_grid1 import make_axes_locatable

color_cycle = cycler('color', ['b', 'g', 'r'])


def plot_figure_16(inputs, outputs, threads, config, params, wildcards):
    # load data
    C_par = inputs["C_par"]
    C_root, ext = path.splitext(inputs['C'][0])
    C = mat.coherencyMatrix(C_root, C_par, gamma=True, bistatic=True,
                            basis='lexicographic').to_monostatic()

    C_zoom, rplot, azplot, (mx_r_zoom, mx_az_zoom), res_dict, r_vec, az_vec = cf.ptarg(C, params['ref']['ridx'],
                                                                                       C.azidx_dec(
                                                                                           params['ref']['azidx']),
                                                                                       rwin=10, azwin=10, sw=(3,3))

    max_VV = cf.dB(np.abs(C_zoom[mx_r_zoom, mx_az_zoom, 2, 2]))
    max_HV = cf.dB(np.abs(C_zoom[mx_r_zoom, mx_az_zoom, 1, 1]))

    r_max = r_vec[mx_r_zoom]
    az_max = az_vec[mx_az_zoom]

    # hv for plotting arrow
    hv_az_plot = (max_HV, az_max)
    hv_r_plot = (r_max, max_HV)
    # vv for plotting arrow
    vv_az_plot = ( max_VV, az_max)
    vv_r_plot = (r_max, max_VV)

    arrow_prop = dict(arrowstyle="<->", )

    plt.style.use(inputs['style'])
    fig_w, fig_h = plt.rcParams['figure.figsize']
    f = plt.figure(figsize=(2 * fig_w, 2 * fig_w))
    ax_pauli = f.add_subplot(111)
    # Axes divider
    divider = make_axes_locatable(ax_pauli)
    ax_r = divider.append_axes('top', size='40%', sharex=ax_pauli, pad=0.1)
    ax_az = divider.append_axes('right', size='40%', sharey=ax_pauli, pad=0.1)
    pauli = vf.pauli_rgb(C_zoom.diagonal(axis1=-2, axis2=-1)[:, :, [2, 1, 0]])
    ext = [r_vec.min(), r_vec.max(), az_vec.min(), az_vec.max(), ]

    # plot pauli rgb
    ax_pauli.imshow(pauli.transpose([1, 0, 2]), extent=ext, aspect=1)
    ax_pauli.xaxis.set_label_text(r'Range [$m$]')
    ax_pauli.yaxis.set_label_text(r'Azimuth [$^\circ$]')
    # plot range response
    ax_r.set_prop_cycle(color_cycle)
    ax_r.plot(r_vec, cf.dB(np.abs(rplot.diagonal(axis1=-2, axis2=-1))))
    ax_r.xaxis.set_visible(False)
    ax_r.yaxis.set_label_text('Intensity [dB]')
    # add arrow for berebbere
    ax_r.annotate("", xy=hv_r_plot, xytext=vv_r_plot, xycoords='data', arrowprops=arrow_prop, textcoords='data')
    lim = [-25, 40]
    ax_r.set_ylim(lim)
    ax_r.grid()
    # plot azimuth response
    ax_az.set_prop_cycle(color_cycle)
    ax_az.yaxis.set_visible(False)
    ax_az.plot(cf.dB(np.abs(azplot.diagonal(axis1=-2, axis2=-1))), az_vec, )
    ax_az.set_xlim(lim)
    ax_az.xaxis.set_label_text('Intensity [dB]')
    ax_az.annotate("", xy=hv_az_plot, xytext=vv_az_plot, xycoords='data', arrowprops=arrow_prop, textcoords='data')

    ax_az.grid()
    plt.show()
    f.savefig(outputs[0])


    # #Geocode
    # LUT = geo.GeocodingTable(inputs.dem_seg_par, inputs.LUT)
    #
    # #Compute decomposition
    # print(C_cal.basis)
    # H, anisotropy, alpha, beta_m, p, w = C_cal.cloude_pottier()
    # #Geocode
    # alpha_geo  = LUT.geocode_data(alpha)
    # H_geo = LUT.geocode_data(H)
    # #RGB composite
    # plot_kw = {'k':0.2, 'sf':1.2, 'coherence':False, 'min_val':0, 'max_val':np.pi/2, 'black_background':False, 'type':'increasing'}
    # span_geo = LUT.geocode_data(C_cal.span())
    # rgb, norm, pal = vf.dismph(span_geo*np.exp(1j*alpha_geo), **plot_kw)
    # pal, pal_ext = vf.dismph_palette(span_geo*np.exp(1j*alpha_geo), **plot_kw)
    # #create figure
    # plt.style.use(inputs['style'])
    # fig_w, fig_h = plt.rcParams['figure.figsize']
    # ext2=LUT.get_geocoded_extent(alpha)
    # ext2 = [ext2[-2], ext2[-1], ext2[0], ext2[1]]
    # #Add axis
    # f1 =  plt.figure(figsize
    #                  =(2 * fig_w,  2 * fig_w))
    # gs = gridspec.GridSpec(*(2, 3), height_ratios=[0.95, 0.05])
    # ax1 = f1.add_subplot(gs[0, ::])
    # cax = f1.add_subplot(gs[1,1])
    # ax1.imshow(rgb.transpose(1, 0, 2), extent=ext2, aspect=vf.fixed_aspect(ext2, 1))
    # ax1.axis('off')
    # cax.imshow(pal, extent=pal_ext, aspect=vf.fixed_aspect(pal_ext, 1))
    # cax.set_ylabel(r'$\alpha$')
    # cax.set_xlabel(r'Intensity')
    # cax.grid(b=False)
    # cax.set_yticks([0, np.pi/4, np.pi/2])
    # cax.set_xticks([0, pal_ext[1]/2, pal_ext[1]])
    # cax.set_yticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$"])
    # gs.update(hspace=0.05)
    # f1.subplots_adjust(left=0, right=1)
    # f1.savefig(outputs['alpha_fig'])


plot_figure_16(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
               snakemake.wildcards)
