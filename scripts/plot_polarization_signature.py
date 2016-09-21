


import pyrat.core.matrices as mat
import pyrat.core.polfun as pf
import pyrat.core.corefun as cf
import matplotlib.pyplot as plt
import pyrat.visualization.visfun as vf
from mpl_toolkits.mplot3d import Axes3D
import pyrat.gpri_utils.calibration as cal
import numpy as np
import csv
def plot_signature(inputs, outputs, threads, config, params, wildcards):
    with plt.style.context(config['style']):
        C = mat.coherencyMatrix(params['C_root'],  params['C_root'] + '.par', gamma=True, bistatic=True, basis='lexicographic').boxcar_filter([3,3])
        #Find the maxium by interpolation
        #Compute the ptarg response for the C matrix
        ptarg_zoom_C, rplot_C, azplot_C, mx_pos, resolution_dict = cf.ptarg(C[:,:,:,:], float(wildcards['ridx']), float(wildcards['azidx']), azwin=20, rwin=10, sw=(2,4))
        ptarg_zoom_C = mat.coherencyMatrix(ptarg_zoom_C, basis='lexicographic', bistatic=True)
        C_tri = ptarg_zoom_C[mx_pos]
        co_sig, x_sig, psi, chi = pf.pol_signature(C_tri.to_monostatic(), n_points=300)
        pp, cc = (np.rad2deg(psi), np.rad2deg(chi))
        co_fig = plt.figure()
        co_ax = co_fig.add_subplot(121,projection='3d')
        # x_ax = f.add_subplot(122,projection='3d',sharex=co_ax, sharey=co_ax)
        co_ax.plot_surface(pp,cc,co_sig, cmap='RdBu_r',lw=0.2)
        co_ax.auto_scale_xyz([0,180],[-45,45], [0,1])
        # co_ax.set_title(r'Copolarized response')
        co_ax.set_ylabel(r' $\chi$ [deg]')
        co_ax.set_xlabel(r' $\psi$ [deg]')
        co_ax.set_zlabel(r'Normalized response')
        x_fig = plt.figure()
        x_ax = x_fig.add_subplot(121,projection='3d', sharex=co_ax, sharey=co_ax)
        x_ax.plot_surface(pp,cc,x_sig,cmap='RdBu_r',lw=0.2)
        x_ax.auto_scale_xyz([0,180],[-45,45], [0,1])
        # x_ax.set_title(r'Cross-polarized response')
        x_ax.set_ylabel(r'$\chi$ [deg]')
        x_ax.set_xlabel(r'$\psi$ [deg]')
        x_ax.set_zlabel(r'Normalized response')
        vf.format_axes(co_ax)
        vf.format_axes(x_ax)
        dist = 15
        xticks = [0,45,90,135,180]
        yticks = [-45,0,45]
        co_ax.dist = dist
        co_ax.set_xticks(xticks)
        co_ax.set_yticks(yticks)
        x_ax.dist = dist
        x_ax.set_xticks(xticks)
        x_ax.set_yticks(yticks)
        for ax in x_ax, co_ax:#ugly hack to change axis width
            ax.zaxis.labelpad = 0.5
            ax.yaxis.labelpad = 0.5
            ax.xaxis.labelpad = 0.5
            ax.w_xaxis.gridlines.set_lw(0.5)
            ax.w_yaxis.gridlines.set_lw(0.5)
            ax.w_zaxis.gridlines.set_lw(0.5)
        co_fig.savefig(outputs['co_sig_plot'],bbox_inches='tight', pad_inches=0)
        x_fig.savefig(outputs['x_sig_plot'],bbox_inches='tight', pad_inches=0)
        plt.tight_layout()

        # ptarg_zoom_VV, rplot_VV, azplot_VV, crap = cf.ptarg(C[:,:,3,3],int(wildcards['ridx']), int(wildcards['azidx']),azwin=100)
        # ptarg_zoom_HV, rplot_HV, azplot_HV, crap = cf.ptarg(C[:, :, 1, 1], int(wildcards['ridx']),int(wildcards['azidx']), azwin=100)
        db = lambda arr: 10 * np.log10(np.abs(arr))
        plot_fig_az = plt.figure()
        plot_ax_az = plot_fig_az.add_subplot(1,1,1)
        plot_ax_az.plot(db(azplot_C[:,0,0]),color='b', label='HH')
        plot_ax_az.plot(db(azplot_C[:,3,3]),color='r', label='VV')
        plot_ax_az.plot(db(azplot_C[:,1,1]), color='g', label='HV')
        # plot_ax_az.set_ylim(-5, 100)
        plot_ax_az.set_xlabel('azimuth samples')
        plot_ax_az.set_ylabel('power [dB]')
        leg = plot_ax_az.legend()
        leg.get_frame().set_linewidth(0.2)
        plot_fig_r = plt.figure()
        plot_ax_r = plot_fig_r.add_subplot(1,1,1)
        plot_ax_r.plot(cf.dB(rplot_C[:,0,0]),color='b', label='HH')
        plot_ax_r.plot(cf.dB(rplot_C[:,3,3]),color='r', label='VV')
        plot_ax_r.plot(cf.dB(rplot_C[:,1,1]), color='g', label='HV')
        plot_ax_r.set_xlabel('range samples')
        plot_ax_r.set_ylabel('power [dB]')
        leg = plot_ax_r.legend()
        leg.get_frame().set_linewidth(0.2)
        plot_fig_az.savefig(outputs['azplot'])
        plot_fig_r.savefig(outputs['rplot'])
        #Now estimate calibration parameters
        r_sl = C.r_vec[int(wildcards['ridx'])]
        HHVV_phase = np.rad2deg(np.angle(C_tri[0,3]))
        f = (C_tri[3,3].real / C_tri[0,0].real)**(1/4.0)
        purity = db(ptarg_zoom_C[mx_pos + (0,0)].real) - db( ptarg_zoom_C[mx_pos + (1,1)].real)
        RCS = db(ptarg_zoom_C[mx_pos + (0,0)].real)
        str = "{HHVV_phase},{f}".format(HHVV_phase=HHVV_phase, f=f)
        with open(params['est_parameters'],'w+') as of:
            tabwrite = csv.writer(of, delimiter=',')
            tabwrite.writerow(['HH-VV phase imbalance', 'HH-VV amplitude imbalance', 'Polarisation purity', 'RCS', 'slant range'])
            tabwrite.writerow([HHVV_phase, f, purity, RCS, r_sl])
plot_signature(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params, snakemake.wildcards)