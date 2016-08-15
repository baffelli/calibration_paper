import os
import glob
import pyrat.fileutils.gpri_files as gpf
import pyrat
import pyrat.gpri_utils
from types import SimpleNamespace
import re
import json






#def reflector_list_from_config():
#	reflectors = config['list_of_reflectors']
#	azidx = [ref[0][1] for ref in reflectors]
#	azidx_dec = [int(ref[0][1]/int(config['range_compression']['dec'])) for ref in reflectors]
#	ridx = [ref[0][0] for ref in reflectors]
#	undec_list = [ (r, az) for r, az in zip(ridx, azidx)]
#	dec_list = [ (r, az) for r, az in zip(ridx, azidx_dec)]
#	return undec_list, dec_list
#Heading
config['geocoding']= {}
config['geocoding']['DEM'] = 'geo/swissALTI3D_2016_Clip.dem'
config['geocoding']['ref_mli_par'] = 'slc/20140321_123707_AAAl.mli.par'
config['geocoding']['scan_heading'] = -70
#Include the necessary files
#include: '../Code/raw_to_slc.snake'
#include: '../Code/slc_to_calibrated_c.snake'
include: '../Code/geocoding.snake'
#Set the working directory to be the current directory
workdir: './'



#Set the rule for geocoding


rule all:
    input:
        "slc/20140321_123707_AAAl.mli_gc.tif",
        "geo/Chutze.inc.tif",
        "geo/Chutze.sim_sar.tif",
        'geo/Chutze.ls_map.tif'


#rule geocoding_table:
#    input:
#        dem = 'geo/swissALTI3D_2016_Clip.dem',
#        dem_par = 'geo/swissALTI3D_2016_Clip.dem_par',
#        reference_mli_par = "slc/{slcname}.par",
#    output:
#        lut = "geo/{slcname}.gpri_to_dem",
#        dem_seg = "geo/{slcname}.dem_seg",
#        dem_seg_par = "geo/{slcname}.dem_seg.par",
#        sim_sar = "geo/{slcname}.sim_sar",
#        lv_theta = "geo/{slcname}.lv_theta",
#        lv_phi = "geo/{slcname}.lv_phi",
#        u = "geo/{slcname}.u",
#        v = "geo/{slcname}.v",
#        inc = "geo/{slcname}.inc",
#        psi = "geo/{slcname}.psi",
#        pix = "geo/{slcname}.pix",
#        ls_map = "geo/{slcname}.ls_map"
#    log:
#        "logs/gc.log"
#    params:
#        lat_ovr = 2,
#        lon_ovr = 2,
#    run:
#        #Determine tiepoints
#        scan_heading = config['geocoding']['scan_heading']
#        cmd = "set_value {{input.reference_mli_par}} {{input.reference_mli_par}} GPRI_scan_heading {scan_heading}".format(scan_heading=scan_heading)
#        shell(cmd)
#        cmd = "gc_GPRI_map {input.reference_mli_par} {input.dem_par} {input.dem} {output.dem_seg_par} {output.dem_seg} {output.lut} {params.lat_ovr} {params.lon_ovr} {output.sim_sar} {output.lv_theta} {output.lv_phi} {output.u} {output.v} {output.inc} {output.psi} {output.pix} {output.ls_map} - > {log}"
#        shell(cmd)
#
##############################################################################
## Geocode using a reference table generated using the methods in the geocoding snakefile
rule geocode_simple:
    input:
        lut = "geo/" + "Chutze.gpri_to_dem",
        dem_seg_par = 'geo/Chutze.dem_seg.par',
        data = "{path_to}/{slcname}",
        reference_mli_par = config['geocoding']['ref_mli_par']
    output:
        geocoded_file = "{path_to}/{slcname}_gc"
    params:
    log:
        "logs/gc.log"
    run:
        dem_width = gpf.get_width(input.dem_seg_par)
        data_width = gpf.get_width(input.reference_mli_par)
        dem_par_dict = gpf.par_to_dict(input.dem_seg_par)
        par_dict = gpf.par_to_dict(input.reference_mli_par)
        nlines = dem_par_dict['nlines']
        filetype = gpf.gamma_datatype_code_from_extension(input.data)
        cmd = "geocode_back  {{input.data}} {data_width} {{input.lut}} {{output.geocoded_file}} {out_width} {nlines} 0 {dtype} >> {{log}}".format(data_width=data_width,
        out_width=dem_width, nlines=nlines, dtype=filetype)
        shell(cmd)


rule mli:
    input:
        slc = "{slcname}.slc",
        slc_par = "{slcname}.slc.par",
    output:
        mli = "{slcname}.mli",
        mli_par = "{slcname}.mli.par",
    shell:
        "multi_look {input.slc} {input.slc_par} {output.mli} {output.mli_par} 2 1"

rule to_bmp:
    input:
        mli = "{slcname}.mli",
        mli_par = "{slcname}.mli.par",
    output:
        bmp = '{slcname}.mli.bmp'
    run:
        mli_par = gpf.par_to_dict(mli_par)
        width = mli_par['width']
        bmp_cmd = "raspwr {input.mli} " + str(width) + " - - - - - - - {output.bmp}"


rule to_gt:
	input:
		file = "{filename}",
		dem_seg_par = "geo/Chutze.dem_seg.par",
	output:
		gt = "{filename}.tif"
	run:
		dt_code = gpf.gt_mapping_from_extension(input.file)
		gt_cmd = "data2geotiff {input.dem_seg_par} {input.file} " + str(dt_code) +  " {output.gt}"
		shell(gt_cmd)


###list of reflectors
#refl_list, refl_list_dec = reflector_list_from_config()
##prefix for all figures
#figures_prefix = 'outputs/img/'  + "2014-2016_GPRI_HIL_Calibration" + '_' + "20140910" + '_' + "144113"
#
##prefix with decimated location
#figures_prefix_with_indices_dec = expand(figures_prefix + "_{ridx}_{azidx}", zip, ridx=[ref[0] for ref in refl_list],
#azidx=[int(ref[1]/config['range_compression']['dec'])  for ref in refl_list])
#
#figures_prefix_with_indices = expand(figures_prefix + "_{ridx}_{azidx}", zip, ridx=[ref[0] for ref in refl_list],
#azidx=[int(ref[1])  for ref in refl_list])
#
##Polarisation signatures
#signatures = expand("{prefixes}_{rx}_{proc}_signature.pdf", proc=['cal', 'flat'], rx='l', prefixes=figures_prefix_with_indices_dec)
#
##Oversampled phase response
#ptarg_mph = expand("{prefixes}_{chan}_{proc}_mph_plot.pdf",proc=['corr', 'desq', 'chan'], chan=['AAAl','BBBl'], prefixes=figures_prefix_with_indices)
#
#phase_plot = expand(figures_prefix + "_{chan}_{processing}_phase_plot.pdf", processing=['coreg', 'corr'], chan=['AAAl', 'BBBl'])
#
#cp_params = expand(figures_prefix + "_l_{type}_gc_{param}.pdf", param=['H', 'alpha', 'pauli'], type=['cal', 'flat'])
#
#
##Making paper
#paper = 'doc/calibration_paper.pdf'
##Combined phase response
#
#rule all:
#	input:
#	    signatures,
#	    ptarg_mph,
#	    phase_plot,
#	    cp_params,
#	    paper
#		HH_VV_phase_plot,
#		sig,
#		squint_plot,
#		ptarg_mph,
#		png_figures,
#		ant_fit,
#		gc_pauli,
#		azimuth_phase_plot,
#		ant_pat,
#		cp_params

##############################################################################
## Moves the figure
#rule move_fig:
#    input:
#        fig = "{location}/{figname}"
#    output:
#        fig = config['figdir'] + {figname}




###############################################################################
### Measure the phase center position and the squint	rate (the two most important
### calibration parameters)
#rule antenna_parameters:
#	input:
#		slc = "{location}/slc_coreg/{date}_{time}_{chan}.slc",
#		ant_function = "{pyrat}measure_phase_center.py".format(pyrat=PYRATDIR),
#		style = config['style']
#	output:
#		fig = "{FIGDIR}/{{location}}_{{date}}_{{time}}_{{ridx}}_{{azidx}}_{{chan}}_phase_fit.{{png}}".format(FIGDIR=FIGDIR)
#	params:
#	    ant_par = "{location}/{date}/ant_par/{date}_{time}_{chan}_{ridx}_{azidx}.ant_par",
#	shell:
#		"""
#		{input.ant_function} {input.slc} {input.slc}.par {wildcards.ridx} {wildcards.azidx}  {params.ant_par} {output.fig} --unwrap
#		"""




###############################################################################
### Produce the phase as an image
#rule dismph:
#	input:
#		ifgram = "{name}",
#		ifgram_par = "{name}.diff_par",
#		function = "{pyrat}mph_image.py".format(pyrat=config['pyrat']),
#	output:
#		ifgram_image = "{name}.mph.{ext}"
#	run:
#		width = int(get_width(input.ifgram_par))
#		cmd = "{{input.function}} {{input.ifgram}} {width} 0 - - 0.1 {{output.ifgram_image}}".format(width=width)
#		shell(cmd)


################################################################################
#### Plot oversampled response
#rule plot_mph:
#	input:
#		slc = "slc_{processing}/{slcname}_{chan}.slc",
#		style = config['style']
#	output:
#		plot = "outputs/img/{slcname}_{ridx}_{azidx}_{chan}_{processing}_mph_plot.{ext}",
#	script:
#		"plotting/plot_mph.py"
#
###############################################################################
### Display polarization signature
#rule pol_signature:
#    input:
#        C = expand("cov_{{type}}/{{slcname}}_{{rx}}.c{i}{j}",i=[0,1,2,3],j=[0,1,2,3]),
#    output:
#        x_sig_plot = "outputs/img/{slcname}_{ridx}_{azidx}_{rx}_{type}_signature_x.{ext}",
#        co_sig_plot = "outputs/img/{slcname}_{ridx}_{azidx}_{rx}_{type}_signature.{ext}",
#        azplot = "outputs/img/{slcname}_{ridx}_{azidx}_{rx}_{type}_azplot.{ext}",
#        rplot = "outputs/img/{slcname}_{ridx}_{azidx}_{rx}_{type}_rplot.{ext}",
#    params:
#        C_root = "cov_{type}/{slcname}_{rx}",
#        est_parameters = "outputs/img/{slcname}_{ridx}_{azidx}_{rx}_{type}_cal_params.csv",
#    script:
#        "plotting/plot_polarization_signature.py"
#
#
###############################################################################
### Plot phase response before and after the correction
#rule plot_phase_response:
#	input:
#		slc = "slc_{processing}/{slcname}_{chan}.slc",
#		style = config['style'],
#	output:
#		plot = "outputs/img/{slcname}_{chan}_{processing}_phase_plot.{ext}"
#	params:
#		ws=50
#	script:
#		"plotting/plot_azimuth_phase.py"
#
###############################################################################
### Geocode using a reference table generated using the methods in the geocoding snakefile
#rule geocode_covariance_with_reference:
#    input:
#        lut = "geo/" + config['gc_reference'],
#        dem_seg_par = "geo/" + config['gc_reference'].replace("gpri_to_dem", "dem_seg.par"),
#        data = "{path_to}/{filename}.{ext}",
#        reference_mli_par = "slc_corr/" + config['gc_reference'].replace(".gpri_to_dem", "_AAAl.slc_dec.par"),
#    output:
#        geocoded_file = "{path_to}/{filename}.{ext}_gc"
#    params:
#    log:
#        "logs/gc.log"
#    run:
#        dem_width = gpf.get_width(input.dem_seg_par)
#        data_width = gpf.get_width(input.reference_mli_par)
#        dem_par_dict = gpf.par_to_dict(input.dem_seg_par)
#        par_dict = gpf.par_to_dict(input.reference_mli_par)
#        nlines = dem_par_dict['nlines']
#        filetype = gpf.gamma_datatype_code_from_extension(input.data)
#        cmd = "geocode_back  {{input.data}} {data_width} {{input.lut}} {{output.geocoded_file}} {out_width} {nlines} 0 {dtype} >> {{log}}".format(data_width=data_width,
#        out_width=dem_width, nlines=nlines, dtype=filetype)
#        shell(cmd)
#
#ruleorder: polcal > plot_polarimetric_analysis > paper
################################################################################
#### Plot Polarimetric analysis
#rule plot_polarimetric_analysis:
#    input:
#        C = expand("cov_{{type}}/{{slcname}}_{{rx}}.c{i}{j}",i=[0,1,2,3],j=[0,1,2,3]),
#        C_gc = expand("cov_{{type}}/{{slcname}}_{{rx}}.c{i}{j}_gc",i=[0,1,2,3],j=[0,1,2,3]),
#        ls_map = "geo/{slcname}.ls_map",
#        lut = "geo/{slcname}.gpri_to_dem",
#        dem_seg_par = "geo/{slcname}.dem_seg.par",
#    params:
#        cov_par =  "cov_{type}/{slcname}_{rx}.par",
#        C_root = "cov_{type}/{slcname}_{rx}"
#    output:
#        alpha = "outputs/img/{slcname}_{rx}_{type}_gc_alpha.{ext}",
#        H = "outputs/img/{slcname}_{rx}_{type}_gc_H.{ext}",
#        pauli = "outputs/img/{slcname}_{rx}_{type}_gc_pauli.{ext}",
#    script:
#        "plotting/plot_polarimetric_analysis.py"
#
#
#
################################################################################
#### Remake the paper when sections or figures change
#rule paper:
#    input:
#        sections = glob.glob('doc/sections/*.tex'),
#        main_paper = 'doc/calibration_paper.tex',
#        figures = glob.glob('outputs/img/*.pdf')
#    output:
#        paper_pdf = 'doc/calibration_paper.pdf'
#    shell:
#        """
#            cd doc
#            pdflatex $(basename {input.main_paper})
#            bibtex $(basename {input.main_paper} .tex)
#            pdflatex $(basename {input.main_paper})
#        """

###############################################################################
### Plot Pauli RGB
#rule plot_gc_pauli:
#    input:
#        C = expand("cov_{{type}}/{{slcname}}_{{rx}}.c{i}{j}",i=[0,1,2,3],j=[0,1,2,3]),
#        C_gc = expand("cov_{{type}}/{{slcname}}_{{rx}}.c{i}{j}_gc",i=[0,1,2,3],j=[0,1,2,3]),
#        ls_map = "geo/{slcname}.ls_map",
#        cov_par = "cov_{type}/{slcname}_{rx}.par",
#        lut = "geo/{slcname}.gpri_to_dem",
#        dem_seg_par = "geo/{slcname}.dem_seg.par",
#    params:
#        C_root = "cov_{type}/{slcname}_{rx}"
#    output:
#        image = "{FIGDIR}/{{location}}_{{date}}_{{time}}_{{rx}}_{{type}}_gc_pauli.{{ext}}".format(FIGDIR=FIGDIR),
#    script:
#        "src/plotting/plot_reflectors.py"
#
#

#
#
#
###############################################################################
### Plot raw data
#rule plot_raw:
#	input:
#		raw = "{location}//desq/{date}_{time}_{chan}.raw",
#		style = config['style']
#	output:
#		image = "{FIGDIR}/{{location}}_{{date}}_{{time}}_{{chan}}_{{azidx}}_raw_plot.{{ext}}".format(FIGDIR=FIGDIR),
#	script:
#		"src/plotting/display_raw_samples.py"
#
#

#
#
###############################################################################
### Plot Phase difference before and after the correction
#rule plot_gc_phase_difference:
#    input:
#        data = "{location}//cov_{proc_type}/{date}_{time}{rx}.c{i}{j}_gct",
#        dem_seg_par = "{location}//geo/{date}_{time}.dem_seg.par",
#        style = config['style']
#    output:
#        image = "{FIGDIR}/{{location}}_{{date,\d+}}_{{time,\d+}}{{rx}}_{{i,\d}}{{j,\d}}_{{proc_type}}_gc_phase.{{ext}}".format(FIGDIR=FIGDIR),
#    script:
#        'src/plotting/plot_phase_difference_geocoded.py'
#
#
#
###############################################################################
### Plot Patterns
#rule plot_pattern:
#	input:
#		HH = expand("/home/baffelli/PhD/trunk/Calibration/src/HH_{freq}GHz.csv",freq=['171','172','173']),
#		VV = expand("/home/baffelli/PhD/trunk/Calibration/src/VV_{freq}GHz.csv",freq=['171','172','173']),
##        chan = "{location}/{date}/desq/{date}_{time}_{chan}.raw"
#	output:
#		plot =  "{FIGDIR}/antenna_pattern.{{ext}}".format(FIGDIR=FIGDIR),
#	script:
#		"src/plotting/plot_patterns.py"
#
#
#
###############################################################################
### Compute squint factor
#rule compute_squint_parameters:
#	input:
#		raw = "{location}//raw_{proc}/{date}_{time}_{chan}.raw",
#		raw_par = "{location}//raw_{proc}/{date}_{time}_{chan}.raw_par",
#        HH = expand("/home/baffelli/PhD/trunk/Calibration/src/HH_{freq}GHz.csv",freq=['171','172','173']),
#        VV = expand("/home/baffelli/PhD/trunk/Calibration/src/VV_{freq}GHz.csv",freq=['171','172','173']),
#		code = "{code_dir}fit_squint.py".format(code_dir=code_dir),
#	output:
#		#squint_vector = "{location}/{date}/{proc}/{date}_{time}_{chan}_{ridx}_{azidx}.squint" ,
#		squint_fig = "{FIGDIR}/{{location}}_{{date}}_{{time}}_{{ridx}}_{{azidx}}_{{chan}}_{{proc}}_squint_plot.{{ext}}".format(FIGDIR=FIGDIR),
#		antenna_pattern_plot =  "{FIGDIR}/{{location}}_{{date}}_{{time}}_{{ridx}}_{{azidx}}_{{chan}}_{{proc}}_squint_plot_antenna_pattern.{{ext}}".format(FIGDIR=FIGDIR)
#	params:
#		plot_flag= lambda wildcards: '--plot_model' if wildcards.proc == 'chan' else '',
#		pattern_name = lambda wildcards: "HH" if wildcards.chan[0:3] == 'AAA' else "VV",
#		squint_par = "{location}//raw_{proc}/{date}_{time}_{ridx}_{azidx}_{chan}.squint_par",
#	script:
#		  "src/fit_squint.py"
#
#
#
##ruleorder: plot_phase_response > pdflatex
#
###############################################################################
### Make a rgb image
#rule rgb:
#	input:
#		r = '{location}//{output_folder}/{date}_{time}_AAA{rx}.mli',
#		g = '{location}//{output_folder}/{date}_{time}_ABB{rx}.mli',
#		b = '{location}//{output_folder}/{date}_{time}_BBB{rx}.mli',
#	output:
#		ras = '{location}//{output_folder}/{date}_{time}_rgb_{rx}.bmp'
#	shell:
#		"""
#		width=$(get_value {input.r}.par range_samples)
#		ras3pwr {input.r} {input.g} {input.b} $width - - - - 2 0 2 1 {output.ras}
#		"""





#Copies the png to the png folder
rule png:
    output:
        "{{name}}.png"
    input:
        "{{name}}.pdf"
    params:
        res = lambda wildcards: 1200 if '_gc' in wildcards.name else 200
    shell:
        "convert -density {params.res} -channel RGBA {input} {output} "






